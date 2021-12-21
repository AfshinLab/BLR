from io import StringIO
from blr.utils import parse_fai, FastaIndexRecord, chromosome_chunks, symlink_relpath, generate_chunks, get_bamtag
from blr.utils import calculate_N50, parse_filters, NaibrSV, parse_naibr_tsv
from pathlib import Path
import os
import pytest

from .test_tagbam import build_read


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


def test_generate_chunks(tmp_path):
    fai = tmp_path / "ref.fasta.fai"
    s = "A\t700\t0\t0\t0\n"   # Phased
    s += "B\t400\t0\t0\t0\n"  # Phased
    s += "C\t300\t0\t0\t0\n"  # Phased
    s += "D\t200\t0\t0\t0\n"  # Primary
    s += "E\t200\t0\t0\t0\n"  # Primary
    s += "F\t100\t0\t0\t0\n"
    with open(fai, "w") as f:
        f.write(s)

    all_contigs = {"A", "B", "C", "D", "E", "F"}
    phasing_contigs = {"A", "B", "C"}
    contigs_skipped = "F"
    primary_contigs = all_contigs - set(contigs_skipped)

    chunks = generate_chunks(str(tmp_path / "ref.fasta"),
                             size=700,
                             phasing_contigs_string=",".join(phasing_contigs),
                             contigs_skipped=contigs_skipped)

    all_chunks_contigs = {c.name for chunk in chunks["all"] for c in chunk}
    phased_chunks_contigs = {c.name for chunk in chunks["phased"] for c in chunk}
    primary_chunks_contigs = {c.name for chunk in chunks["primary"] for c in chunk}

    assert all_chunks_contigs == all_contigs
    assert primary_chunks_contigs == primary_contigs
    assert phased_chunks_contigs == phasing_contigs

    assert primary_chunks_contigs.issubset(all_chunks_contigs)
    assert phased_chunks_contigs.issubset(primary_chunks_contigs)


def test_symlink_relpath(tmpdir):
    path_a = Path(tmpdir / "folder" / "fileA.txt")
    path_b = Path(tmpdir / "folder" / "fileB.txt")
    path_a.parent.mkdir()
    path_a.touch()
    symlink_relpath(path_a, path_b)
    assert path_b.exists()
    assert path_b.is_symlink()
    assert path_b.samefile(path_a)
    assert os.readlink(path_b) == path_a.name


def test_get_bamtag():
    barcode = "ATGCATGC"
    read1 = build_read(name="read1", barcode=barcode)
    read1_barcode = get_bamtag(read1, "BX")
    assert read1_barcode == barcode

    read2 = build_read(name="read2")
    read2_barcode = get_bamtag(read2, "BX")
    assert read2_barcode is None

    default_barcode = "DEFAULT"
    read2_barcode = get_bamtag(read2, "BX", default="DEFAULT")
    assert read2_barcode == default_barcode


def test_calculate_N50():
    # Used example from https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics#N50
    # Sum of values is 54, half is 27. 10 + 9 + 8 = 27 --> N50 is 8
    values = [2, 3, 4, 5, 6, 7, 8, 9, 10]
    assert calculate_N50(values) == 8


def test_parse_filters():
    filter_command_ref = " -filter 'QUAL < 15' --filter-name 'lowQUAL'" \
                         " -filter 'MQRankSum > 6.0' --filter-name 'highMQRankSum'" \
                         " -filter 'AF < 0.15' --filter-name 'snpLowAF'"
    filter_string = "'QUAL < 15','lowQUAL';'MQRankSum > 6.0','highMQRankSum';'AF < 0.15','snpLowAF'"
    filter_command = parse_filters(filter_string)

    assert filter_command == filter_command_ref


def test_naibr_sv_from_string():
    string = "22\t36746940\t22\t36751965\t40.0\t0.0\t+-\t3,3\t258.872434847\tPASS"
    sv = NaibrSV.from_string(string)
    assert sv.pass_filter == "PASS"
    assert sv.break2 == 36_751_965


def test_parse_naibr_tsv():
    s = "Chr1\tBreak1\tChr2\tBreak2\tSplit molecules\tDiscordant reads\tOrientation\tHaplotype\tScore\tPass filter\n"
    s += "22\t36746940\t22\t36751965\t40.0\t0.0\t+-\t3,3\t258.872434847\tPASS\n"
    s += "22\t65743624\t22\t65783121\t33.0\t3.0\t--\t2,2\t41.241421411\tFAIL\n"
    with StringIO(s) as f:
        parser = parse_naibr_tsv(f)
        first = next(parser)
        assert first.pass_filter == "PASS"
        second = next(parser)
        assert second.pass_filter == "FAIL"

        with pytest.raises(StopIteration):
            next(parser)
