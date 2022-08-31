from collections import Counter
import pysam

from blr.cli.tagbam import is_sequence, mode_ema, mode_samtags_underline_separation, get_mode, mode_void


def build_read(name, barcode=None):
    # Read building modified from pysam tests in:
    # https://github.com/pysam-developers/pysam/blob/531cd8569ff5a3ce534956d7795300d52bab4177/tests/AlignedSegment_test.py#L18
    header = pysam.AlignmentHeader.from_references(["chr1"], [10000000])
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.query_sequence = "ATGC" * 10
    a.flag = 0
    a.reference_id = 0
    a.reference_start = 20
    a.mapping_quality = 20
    a.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
    a.next_reference_id = 0
    a.next_reference_start = 200
    a.template_length = 167
    a.query_qualities = pysam.qualitystring_to_array("1234") * 10
    if barcode:
        a.tags = (("BX", barcode),)
    return a


def test_is_sequence():
    assert is_sequence("ATGCTAGTC")
    assert not is_sequence("12318321")


# Typical for bowtie2, bwa, minimap2 mappers.
def test_default_formatted_read():
    name = "myread:11213:21314"
    barcode = "ACTGACTGACTGACTGACTG"
    sequence = "ACTGACTGTCTGACTGCCTG"
    read = build_read(f"{name}_BX:Z:{barcode}_RX:Z:{sequence}")

    modified_read = read.__copy__()
    sample_nr = 3
    mode_samtags_underline_separation(modified_read, 3, "BX", Counter())
    assert modified_read.get_tag("BX") == f"{barcode}-{sample_nr}"
    assert modified_read.get_tag("RX") == sequence
    assert modified_read.query_name == name


def test_default_formatted_read_no_barcode():
    name = "myread:11213:21314"
    read = build_read(name)

    modified_read = read.__copy__()
    mode_samtags_underline_separation(modified_read, 1, "BX", Counter())
    assert not modified_read.has_tag("BX")
    assert not modified_read.has_tag("RX")
    assert modified_read.query_name == name


def test_ema_formatted_read():
    name = "myread:11213:21314"
    barcode = "ACTGACTGACTGACTGACTG"
    barcode_in_tag = barcode[:16] + "-1"
    read = build_read(f"{name}:{barcode}", barcode=barcode_in_tag)

    modified_read = read.__copy__()
    sample_nr = 3
    mode_ema(modified_read, 3, "BX", Counter())
    assert modified_read.get_tag("BX") == f"{barcode}-{sample_nr}"
    assert modified_read.query_name == name


def test_ema_formatted_read_no_barcode():
    name = "myread:11213:21314"
    read = build_read(name)

    modified_read = read.__copy__()
    mode_ema(modified_read, 3, "BX", Counter())
    assert not modified_read.has_tag("BX")
    assert modified_read.query_name == name


def test_get_mode_ema():
    name = "myread:11213:21314"
    barcode = "ACTGACTGACTGACTGACTG"
    barcode_in_tag = barcode[:16] + "-1"
    read = build_read(f"{name}:{barcode}", barcode=barcode_in_tag)

    mode, _ = get_mode(iter([read]), "BX")
    assert mode.__name__ == mode_ema.__name__


def test_get_mode_samtags_underline_separation():
    name = "myread:11213:21314"
    barcode = "ACTGACTGACTGACTGACTG"
    sequence = "ACTGACTGTCTGACTGCCTG"
    read = build_read(f"{name}_BX:Z:{barcode}_RX:Z:{sequence}")

    mode, _ = get_mode(iter([read]), "BX")
    assert mode.__name__ == mode_samtags_underline_separation.__name__


def test_get_mode_void():
    name = "myread:11213:21314"
    read = build_read(name)

    mode, _ = get_mode(iter([read]), "BX")
    assert mode.__name__ == mode_void.__name__
