from pathlib import Path
import shutil
import pysam
import pytest
import dnaio
from xopen import xopen

from blr.__main__ import main as blr_main
from blr.cli.init import init, init_from_dir
from blr.cli.run import run
from blr.cli.config import change_config
from blr.utils import get_bamtag

TESTDATA = Path("blr-testdata")

TESTDATA_BLR_READ1 = TESTDATA / "blr_reads.1.fastq.gz"
TESTDATA_BLR_READ2 = TESTDATA / "blr_reads.2.fastq.gz"
TESTDATA_TENX_READ1 = TESTDATA / "tenx_reads.1.fastq.gz"
TESTDATA_TENX_READ2 = TESTDATA / "tenx_reads.2.fastq.gz"
TESTDATA_TENX_BARCODES = str((TESTDATA / "tenx_barcode_whitelist.txt").absolute())
TESTDATA_STLFR_READ1 = TESTDATA / "stlfr_reads.1.fastq.gz"
TESTDATA_STLFR_READ2 = TESTDATA / "stlfr_reads.2.fastq.gz"
TESTDATA_STLFR_BARCODES = str((TESTDATA / "stlfr_barcodes.txt").absolute())
TESTDATA_TELLSEQ_READ1 = TESTDATA / "tellseq_reads.1.fastq.gz"
TESTDATA_TELLSEQ_READ2 = TESTDATA / "tellseq_reads.2.fastq.gz"
TESTDATA_TELLSEQ_INDEX = str((TESTDATA / "tellseq_index.fastq.gz").absolute())
DEFAULT_CONFIG = "blr.yaml"
REFERENCE_GENOME = str((TESTDATA / "ref.fasta").absolute())
REFERENCE_VARIANTS = str((TESTDATA / "HG002_GRCh38_GIAB_highconf.vcf").absolute())
DB_SNP = str((TESTDATA / "dbSNP.vcf.gz").absolute())


def count_bam_alignments(path):
    with pysam.AlignmentFile(path) as af:
        n = 0
        for _ in af:
            n += 1
    return n


def count_fastq_reads(path):
    with dnaio.open(path) as f:
        n = 0
        for _ in f:
            n += 1
    return n


def count_lariat_fastq_reads(path):
    with xopen(path) as f:
        return sum(1 for line in f if line.startswith("@"))


def count_bam_tags(path, tag):
    with pysam.AlignmentFile(path, 'rb') as af:
        n = 0
        for r in af.fetch(until_eof=True):
            if get_bamtag(r, tag):
                n += 1
    return n


def bam_has_tag(path, tag):
    with pysam.AlignmentFile(path) as file:
        for alignment in file:
            if alignment.has_tag(tag):
                return True

    return False


@pytest.fixture(scope="module")
def _workdir(tmp_path_factory):
    """
    This runs the pipeline using default parameters up to the creation of
    the BAM file ready for variant calling
    """
    path = tmp_path_factory.mktemp(basename="analysis-") / "analysis"
    init(path, TESTDATA_BLR_READ1, "blr")
    change_config(
        path / DEFAULT_CONFIG, [
            ("genome_reference", REFERENCE_GENOME),
            ("chunk_size", "50000")
        ]
    )
    # chromosomes B, C and D end up in the same chunk
    run(workdir=path, targets=[f"chunks/chr{c}.calling.bam.bai" for c in "AB"])
    return path


@pytest.fixture
def workdir(_workdir, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir, path)
    return path


def test_init(tmp_path):
    init(tmp_path / "analysis", TESTDATA_BLR_READ1, "blr")


def test_config(tmp_path):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_BLR_READ1, "blr")
    change_config(workdir / "blr.yaml", [("read_mapper", "bwa")])


def test_default_read_mapper(workdir):
    n_input_fastq_reads = 2 * count_fastq_reads(workdir / "trimmed_barcoded.1.fastq.gz")
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "initialmapping.bam")


# The read mapper will partly determine the output format so we test for different mappers here. Bowtie2, bwa and
# minimap2 all use the same format so only bowtie2 is tested.
@pytest.mark.parametrize("read_mapper", ["bowtie2", "ema"])
def test_trim_blr(workdir, read_mapper):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", read_mapper)]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed, force_run=["trim"])
    assert count_fastq_reads(workdir / trimmed[0]) / count_fastq_reads(TESTDATA_BLR_READ1) > 0.9
    assert count_fastq_reads(workdir / trimmed[1]) / count_fastq_reads(TESTDATA_BLR_READ2) > 0.9


@pytest.mark.skipif(shutil.which("lariat") is None, reason="Lariat not installed")
def test_trim_blr_lariat(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", "lariat")]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed, force_run=["trim"])
    assert count_lariat_fastq_reads(workdir / trimmed[0]) / count_fastq_reads(TESTDATA_BLR_READ1) > 0.9


@pytest.mark.parametrize("read_mapper", ["bowtie2", "ema"])
def test_trim_tenx(tmp_path, read_mapper):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_TENX_READ1, "10x")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("barcode_whitelist", TESTDATA_TENX_BARCODES),
         ("read_mapper", read_mapper)]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed)
    for raw, trimmed in zip((TESTDATA_TENX_READ1, TESTDATA_TENX_READ2), trimmed):
        assert count_fastq_reads(workdir / trimmed) / count_fastq_reads(raw) > 0.9  # More than 90% kept


@pytest.mark.parametrize("read_mapper", ["bowtie2", "ema"])
def test_trim_stlfr(tmp_path, read_mapper):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_STLFR_READ1, "stlfr")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("stlfr_barcodes", TESTDATA_STLFR_BARCODES),
         ("read_mapper", read_mapper)]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed)
    for raw, trimmed in zip((TESTDATA_STLFR_READ1, TESTDATA_STLFR_READ2), trimmed):
        assert count_fastq_reads(workdir / trimmed) / count_fastq_reads(raw) > 0.7


@pytest.mark.parametrize("read_mapper", ["bowtie2", "ema"])
def test_trim_tellseq(tmp_path, read_mapper):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_TELLSEQ_READ1, "tellseq")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("tellseq_index", TESTDATA_TELLSEQ_INDEX),
         ("read_mapper", read_mapper)]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed)
    for raw, trimmed in zip((TESTDATA_TELLSEQ_READ1, TESTDATA_TELLSEQ_READ2), trimmed):
        assert count_fastq_reads(raw) / count_fastq_reads(workdir / trimmed) > 0.9


@pytest.mark.skipif(shutil.which("lariat") is None, reason="Lariat not installed")
def test_trim_tellseq_lariat(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("tellseq_index", TESTDATA_TELLSEQ_INDEX),
         ("read_mapper", "lariat")]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed, force_run=["trim"])
    assert count_lariat_fastq_reads(workdir / trimmed[0]) / count_fastq_reads(TESTDATA_BLR_READ1) > 0.9


non_default_mappers = ["bwa", "minimap2", "ema", "lariat"] if shutil.which("lariat") else ["bwa", "minimap2", "ema"]


@pytest.mark.parametrize("read_mapper", non_default_mappers)
def test_nondefault_read_mappers(tmp_path, read_mapper):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_BLR_READ1, "blr")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME), ("read_mapper", read_mapper)]
    )
    run(workdir=workdir, targets=["initialmapping.bam"])
    n_input_fastq_reads = 2 * count_fastq_reads(workdir / "trimmed_barcoded.1.fastq.gz")
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "initialmapping.bam")


def test_final_compressed_reads_exist(workdir):
    targets = ("reads.1.final.fastq.gz", "reads.2.final.fastq.gz")
    run(workdir=workdir, targets=targets)
    for filename in targets:
        assert workdir.joinpath(filename).exists()


def test_BQSR(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("dbSNP", DB_SNP), ("BQSR", "true"), ("reference_variants", "null"),
         ("variant_caller", "gatk")]
    )
    # Since the workdir fixture creates the calling.bam files already, we need to
    # ensure they are re-created with BQSR applied
    for calling_bam in workdir.joinpath("chunks").glob("*.calling.bam"):
        calling_bam.unlink()
    target = "final.bam"
    run(workdir=workdir, targets=[target])
    with pysam.AlignmentFile(workdir / target) as af:
        # Ensure that ApplyBQSR was run on the file by inspecting the @PG lines in the header
        bqsr_header_entries = [entry for entry in af.header["PG"] if entry["ID"] == "GATK ApplyBQSR"]
        assert bqsr_header_entries


@pytest.mark.parametrize("variant_caller", ["freebayes", "bcftools", "gatk"])
def test_call_variants(workdir, variant_caller):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null"), ("variant_caller", variant_caller)]
    )
    target = "called.vcf"
    run(workdir=workdir, targets=[target])
    assert workdir.joinpath(target).is_file()


def test_filter_variants(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null"), ("filter_variants", "true")]
    )
    target = "called.filtered.vcf"
    run(workdir=workdir, targets=[target])
    assert workdir.joinpath(target).is_file()


def test_plot_figures(workdir):
    target = "figures"
    run(workdir=workdir, targets=[target])
    assert workdir.joinpath(target).is_dir()
    # Check that folder contains 7 PNG images for multiqc.
    assert sum(file.name.endswith("_mqc.png") for file in workdir.joinpath(target).iterdir()) == 7


@pytest.mark.parametrize("haplotype_tool", ["blr", "whatshap"])
def test_haplotag(workdir, haplotype_tool):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null")]
    )
    target = "final.phased.bam"
    run(workdir=workdir, targets=[target])
    assert bam_has_tag(workdir / target, "HP")
    assert bam_has_tag(workdir / target, "PS")
    assert count_bam_tags(workdir / target, "PS") == count_bam_tags(workdir / target, "HP")


def test_version_exit_code_zero():
    with pytest.raises(SystemExit) as e:
        blr_main(["--version"])
    assert e.value.code == 0


def test_init_from_workdir(tmp_path, workdir):
    old_workdir = workdir
    new_workdir = tmp_path / "from_old"

    # Generate all require files is old workdir
    run(workdir=old_workdir, targets=["final.bam", "final.molecule_stats.filtered.tsv"])

    # Initialize new dir based on old and run setup.
    init_from_dir(new_workdir, [old_workdir], "blr")
    change_config(
        new_workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME),
         ("chunk_size", "50000")]
        )
    run(workdir=new_workdir, snakefile="run_anew.smk")
