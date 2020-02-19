import gzip
import pyfaidx


class FastaNotIndexedError(Exception):
    pass


def detect_file_format(path):
    """
    Detect file format and return 'BAM', 'CRAM', 'VCF' or None. None indicates an
    unrecognized file format.

    'VCF' is returned for both uncompressed and compressed VCFs (.vcf and .vcf.gz).
    """
    with open(path, "rb") as f:
        first_bytes = f.read(16)
        if first_bytes.startswith(b"CRAM"):
            return "CRAM"
        if first_bytes.startswith(b"##fileformat=VCF"):
            return "VCF"

    gzip_header = b"\037\213"
    if first_bytes.startswith(gzip_header):
        with gzip.GzipFile(path, "rb") as f:
            first_bytes = f.read(16)
            if first_bytes.startswith(b"BAM\1"):
                return "BAM"
            elif first_bytes.startswith(b"##fileformat=VCF"):
                return "VCF"

    return None


def IndexedFasta(path):
    try:
        f = pyfaidx.Fasta(path, as_raw=True, sequence_always_upper=True, build_index=False)
    except pyfaidx.IndexNotFoundError:
        raise FastaNotIndexedError(path)
    return f


def plural_s(n: int) -> str:
    return "" if n == 1 else "s"
