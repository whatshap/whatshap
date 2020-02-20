import sys
import resource
import logging

from whatshap.bam import (
    AlignmentFileNotIndexedError,
    EmptyAlignmentFileError,
    SampleNotFoundError,
    ReferenceNotFoundError,
)
from whatshap.variants import ReadSetReader, ReadSetError
from whatshap.utils import IndexedFasta, FastaNotIndexedError
from whatshap.core import ReadSet


logger = logging.getLogger(__name__)


class CommandLineError(Exception):
    """An anticipated command-line error occurred. This ends up as a user-visible error message"""

    pass


def open_readset_reader(*args, **kwargs):
    try:
        readset_reader = ReadSetReader(*args, **kwargs)
    except OSError as e:
        raise CommandLineError(e)
    except AlignmentFileNotIndexedError as e:
        raise CommandLineError(
            "The file {!r} is not indexed. Please create the appropriate BAM/CRAM "
            'index with "samtools index"'.format(e)
        )
    except EmptyAlignmentFileError as e:
        raise CommandLineError(
            "No reads could be retrieved from {!r}. If this is a CRAM file, possibly the "
            "reference could not be found. Try to use --reference=... or check you "
            "$REF_PATH/$REF_CACHE settings".format(e)
        )
    return readset_reader


class PhasedInputReader:
    def __init__(
        self,
        bam_paths,
        vcf_readers,
        reference,
        numeric_sample_ids,
        ignore_read_groups,
        **kwargs,  # passed to ReadSetReader constructor
    ):
        # TODO exit stack!
        self._numeric_sample_ids = numeric_sample_ids
        self._fasta = self._open_reference(reference) if reference else None
        self._vcf_readers = vcf_readers
        self._ignore_read_groups = ignore_read_groups

        self._readset_reader = open_readset_reader(
            bam_paths, reference, numeric_sample_ids, **kwargs,
        )

    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self._fasta is not None:
            self._fasta.close()

    @staticmethod
    def _open_reference(path):
        try:
            indexed_fasta = IndexedFasta(path)
        except (IOError, OSError) as e:
            raise CommandLineError("Error while opening FASTA reference file: {}".format(e))
        except FastaNotIndexedError as e:
            raise CommandLineError(
                "An index file (.fai) for the reference FASTA {!r} "
                "could not be found. Please create one with "
                '"samtools faidx".'.format(e)
            )
        return indexed_fasta

    def read(self, chromosome, variants, sample, *, read_vcf=True, regions=None):
        """
        Return a pair (readset, vcf_source_ids) where readset is a sorted ReadSet.

        Set read_vcf to False to not read phased blocks from the VCFs
        """
        readset_reader = self._readset_reader
        for_sample = "for sample {!r} ".format(sample) if not self._ignore_read_groups else ""
        logger.info("Reading alignments %sand detecting alleles ...", for_sample)
        try:
            reference = self._fasta[chromosome] if self._fasta else None
        except KeyError:
            raise CommandLineError(
                "Chromosome {!r} present in VCF file, but not in the reference FASTA {!r}".format(
                    chromosome, self._fasta.filename
                )
            )

        bam_sample = None if self._ignore_read_groups else sample
        try:
            readset = readset_reader.read(chromosome, variants, bam_sample, reference, regions)
        except SampleNotFoundError:
            logger.warning("Sample %r not found in any BAM/CRAM file.", bam_sample)
            readset = ReadSet()
        except ReadSetError as e:
            raise CommandLineError(e)
        except ReferenceNotFoundError:
            if chromosome.startswith("chr"):
                alternative = chromosome[3:]
            else:
                alternative = "chr" + chromosome
            message = "The chromosome {!r} was not found in the BAM/CRAM file.".format(chromosome)
            if readset_reader.has_reference(alternative):
                message += " Found {!r} instead".format(alternative)
            raise CommandLineError(message)

        vcf_source_ids = set()
        if read_vcf:
            # Add phasing information from VCF files, if present
            sample_id = self._numeric_sample_ids[sample]
            for i, vcf_reader in enumerate(self._vcf_readers):
                if chromosome in vcf_reader:
                    variant_table = vcf_reader[chromosome]
                    source_id = readset_reader.n_paths + i
                    vcf_source_ids.add(source_id)
                    for read in variant_table.phased_blocks_as_reads(
                        sample, variants, source_id, sample_id
                    ):
                        readset.add(read)

        # TODO is this necessary?
        for read in readset:
            read.sort()
        readset.sort()

        logger.info(
            "Found %d reads covering %d variants", len(readset), len(readset.get_positions()),
        )
        return readset, vcf_source_ids


def log_memory_usage():
    if sys.platform == "linux":
        memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum memory usage: %.3f GB", memory_kb / 1e6)
