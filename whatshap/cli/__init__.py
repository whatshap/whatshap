from whatshap.bam import AlignmentFileNotIndexedError, EmptyAlignmentFileError
from whatshap.variants import ReadSetReader
from whatshap.utils import IndexedFasta, FastaNotIndexedError


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


def open_reference(path):
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
