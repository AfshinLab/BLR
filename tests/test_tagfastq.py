from collections import Counter
from pytest import raises

from blr.cli.tagfastq import match_template, IUPAC, scramble, map_corrected_barcodes, BarcodeReader

from .utils import tempinput


def test_match_template():
    template = [set(IUPAC[base]) for base in "AWVN"]
    assert match_template("ATCG", template)
    assert not match_template("GTCG", template)
    assert not match_template("ATC", template)


def test_scramble():
    sequences = [
        "AAAAAAAAAAAAAAAAAAAA",
        "AAAAAAAAAAAAAAAATTTT",
        "TTTTTTTTTTTTTTTTTTTT",
        "CCCCCCCCCCCCCCCCTTTT",
    ]
    scrambled_sequences = sequences.copy()
    scramble(scrambled_sequences)
    assert len(scrambled_sequences) == len(sequences)
    assert all(s in scrambled_sequences for s in sequences)
    assert all(s1[:16] != s2[:16] for s1, s2 in zip(scrambled_sequences[:-1], scrambled_sequences[1:]))
    with raises(AssertionError):
        assert all(s1[:16] != s2[:16] for s1, s2 in zip(sequences[:-1], sequences[1:]))


def test_map_barcodes():
    filetext = b"AAAA\t3\tAAAA,AAAT\n" \
               b"TTAA\t2\tTTAG\n" \
               b"CGTG\t1\tCGTG\n"

    with tempinput(filetext) as file:
        corrected_barcodes, _ = map_corrected_barcodes(file, summary=Counter(), mapper="bowtie2",
                                                       template=None)
        assert corrected_barcodes["AAAT"] == "AAAA"
        assert len(corrected_barcodes) == 4
        assert len(set(corrected_barcodes.values())) == 3

        corrected_barcodes, _ = map_corrected_barcodes(file, summary=Counter(), mapper="bowtie2",
                                                       template=None, min_count=2)
        assert len(corrected_barcodes) == 3
        assert len(set(corrected_barcodes.values())) == 2


def test_barcode_parsing_split_header():
    fasta = b">MYHEADER EXTRA\nGATTACCA"
    with tempinput(fasta) as f:
        bcreader = BarcodeReader(f)
        barcode = bcreader.get_barcode("MYHEADER")
        assert barcode == "GATTACCA"


def test_barcode_parsing_single_header():
    fasta = b">MYHEADER\nGATTACCA"
    with tempinput(fasta) as f:
        bcreader = BarcodeReader(f)
        barcode = bcreader.get_barcode("MYHEADER")
        assert barcode == "GATTACCA"
