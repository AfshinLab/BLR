from collections import Counter

from blr.cli.tagfastq import match_template, IUPAC, scramble, parse_corrected_barcodes

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
    for s1, s2 in zip(scrambled_sequences[:-1], scrambled_sequences[1:]):
        assert s1[:16] != s2[:16]


def test_parse_barcodes():
    filetext = b"AAAA\t3\tAAAA,AAAT\n" \
               b"TTAA\t2\tTTAG\n" \
               b"CGTG\t1\tCGTG\n"

    with tempinput(filetext) as file:
        with open(file) as f:
            corrected_barcodes, _ = parse_corrected_barcodes(open_file=f, summary=Counter(), mapper="bowtie2",
                                                             template=None)
            assert corrected_barcodes["AAAT"] == "AAAA"
            assert len(corrected_barcodes) == 4
            assert len(set(corrected_barcodes.values())) == 3

        with open(file) as f:
            corrected_barcodes, _ = parse_corrected_barcodes(open_file=f, summary=Counter(), mapper="bowtie2",
                                                             template=None, skip_singles=True)
            assert len(corrected_barcodes) == 3
            assert len(set(corrected_barcodes.values())) == 2
