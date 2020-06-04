from pathlib import Path
import shutil
import pysam
import pytest
import dnaio

from blr.cli.init import init
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
        path / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME)]
    )
    run(workdir=path, targets=["mapped.calling.bam.bai"])
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
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "mapped.sorted.bam")


def test_trim_blr(workdir):
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed)
    assert count_fastq_reads(trimmed[0]) <= count_fastq_reads(TESTDATA_BLR_READ1)
    assert count_fastq_reads(trimmed[1]) <= count_fastq_reads(TESTDATA_BLR_READ2)


def test_trim_tenx(tmp_path):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_TENX_READ1, "10x")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("barcode_whitelist", TESTDATA_TENX_BARCODES)]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed)
    for raw, trimmed in zip((TESTDATA_TENX_READ1, TESTDATA_TENX_READ2), trimmed):
        assert count_fastq_reads(raw) == count_fastq_reads(workdir / trimmed)


def test_trim_stlfr(tmp_path):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_STLFR_READ1, "stlfr")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("stlfr_barcodes", TESTDATA_STLFR_BARCODES)]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, targets=trimmed)
    for raw, trimmed in zip((TESTDATA_STLFR_READ1, TESTDATA_STLFR_READ2), trimmed):
        assert count_fastq_reads(raw) >= count_fastq_reads(workdir / trimmed)


@pytest.mark.parametrize("read_mapper", ["bwa", "minimap2", "ema"])
def test_nondefault_read_mappers(tmp_path, read_mapper):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_BLR_READ1, "blr")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME), ("read_mapper", read_mapper)]
    )
    run(workdir=workdir, targets=["mapped.sorted.bam"])
    n_input_fastq_reads = 2 * count_fastq_reads(workdir / "trimmed_barcoded.1.fastq.gz")
    assert n_input_fastq_reads <= count_bam_alignments(workdir / "mapped.sorted.bam")


def test_final_compressed_reads_exist(workdir):
    targets = ("reads.1.final.fastq.gz", "reads.2.final.fastq.gz")
    run(workdir=workdir, targets=targets)
    for filename in targets:
        assert workdir.joinpath(filename).exists()


def test_link_reference_variants(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", REFERENCE_VARIANTS)]
    )
    target = "mapped.phaseinput.vcf"
    run(workdir=workdir, targets=[target])
    assert workdir.joinpath(target).is_symlink()


def test_BQSR(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("dbSNP", DB_SNP), ("BQSR", "true"), ("reference_variants", "null"),
         ("variant_caller", "gatk")]
    )
    target = "mapped.sorted.tag.bcmerge.mkdup.mol.filt.BQSR.bam"
    run(workdir=workdir, targets=[target])
    assert workdir.joinpath(target).is_file()


@pytest.mark.parametrize("variant_caller", ["freebayes", "bcftools", "gatk"])
def test_call_variants(workdir, variant_caller):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null"), ("variant_caller", variant_caller)]
    )
    target = "mapped.variants.called.vcf"
    run(workdir=workdir, targets=[target])
    assert workdir.joinpath(target).is_file()


def test_plot_figures(workdir):
    target = "figures/mapped"
    run(workdir=workdir, targets=[target])
    assert workdir.joinpath(target).is_dir()


@pytest.mark.parametrize("haplotype_tool", ["blr", "whatshap"])
def test_haplotag(workdir, haplotype_tool):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null")]
    )
    target = "mapped.calling.phased.bam"
    run(workdir=workdir, targets=[target])
    assert bam_has_tag(workdir / target, "HP")
    assert bam_has_tag(workdir / target, "PS")
    assert count_bam_tags(workdir / target, "PS") == count_bam_tags(workdir / target, "HP")
