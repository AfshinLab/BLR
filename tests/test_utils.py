from io import StringIO
from blr.utils import parse_fai, FastaIndexRecord, chromosome_chunks


def test_parse_fai():
    s = "chr1\t112233\t112\t70\t71\n"
    s += "\n"
    s += "chr2\t445566\t252513167\t70\t71"
    with StringIO(s) as f:
        records = parse_fai(f)
    assert len(records) == 2
    assert records[0] == FastaIndexRecord("chr1", 112233)
    assert records[1] == FastaIndexRecord("chr2", 445566)


def test_chromosome_chunks():
    a = FastaIndexRecord("A", 500)
    b = FastaIndexRecord("B", 500)
    c = FastaIndexRecord("C", 10)
    d = FastaIndexRecord("D", 20)
    e = FastaIndexRecord("E", 150)
    f = FastaIndexRecord("F", 50)
    records = [a, b, c, d, e, f]
    chunks = list(chromosome_chunks(records, size=100))
    assert chunks == [[a], [b], [c, d], [e], [f]]


def test_chromosome_chunks_exact():
    a = FastaIndexRecord("A", 500)
    b = FastaIndexRecord("B", 500)
    chunks = list(chromosome_chunks([a, b], size=1000))
    assert chunks == [[a, b]]
