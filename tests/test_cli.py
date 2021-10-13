from pathlib import Path
import shutil
import pysam
import pytest
import dnaio
from xopen import xopen

from blr.__main__ import main as blr_main
from blr.cli.init import CONFIGURATION_FILE_NAME, init, init_from_dir
from blr.cli.run import run
from blr.cli.config import change_config
from blr.utils import get_bamtag

TESTDATA = Path("blr-testdata")

TESTDATA_DBS_READ1 = TESTDATA / "dbs_reads.1.fastq.gz"
TESTDATA_DBS_READ2 = TESTDATA / "dbs_reads.2.fastq.gz"
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
REFERENCE_VARIANTS = str((TESTDATA / "HG002_GRCh38_GIAB_highconf.vcf.gz").absolute())
DB_SNP = str((TESTDATA / "dbSNP.vcf.gz").absolute())

DEFAULT_SMK_ARGS = ["--notemp", "--show-failed-logs"]


def count_bam_alignments(path):
    with pysam.AlignmentFile(path) as af:
        n = 0
        for a in af:
            if not a.is_secondary and not a.is_supplementary:
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


def chromsome_phased_in_vcf(path, chromosome):
    with pysam.VariantFile(path) as f:
        for rec in f.fetch(chromosome):
            if "PS" in rec.samples[0] and rec.samples[0]["PS"] is not None:
                return True
    return False


def fastq_ema_has_barcodes_grouped(fastq_bin):
    barcodes = [None]
    with xopen(fastq_bin) as f:
        for line in f:
            if line.startswith("@"):
                barcode = line.split(" ")[-1]
                if barcode.startswith("BX") and barcodes[-1] != barcode:
                    barcodes.append(barcode)
    return len(barcodes) == len(set(barcodes))


def fastq_lariat_has_barcodes_grouped(fastq_lariat):
    barcodes = [None]
    with xopen(fastq_lariat) as f:
        entry_row = 9
        for line in f:
            if line.startswith("@") and entry_row == 9:
                entry_row = 0
                continue

            entry_row += 1
            if entry_row == 5 and barcodes[-1] != line.strip():
                barcodes.append(line.strip())
    return len(barcodes) == len(set(barcodes))


@pytest.fixture(scope="module")
def _workdir(tmp_path_factory):
    """
    This runs the pipeline using default parameters up to the creation of
    the BAM file ready for variant calling
    """
    path = tmp_path_factory.mktemp(basename="analysis-") / "analysis"
    init(path, TESTDATA_DBS_READ1, "dbs")
    change_config(
        path / DEFAULT_CONFIG, [
            ("genome_reference", REFERENCE_GENOME),
            ("chunk_size", "50000"),
            ("phasing_contigs", "null"),
            ("heap_space", "1"),
            ("min_mapq", "0"),
            ("fastq_bins", "5"),
        ]
    )
    # chromosomes B, C and D end up in the same chunk
    targets = [f"chunks/chr{c}.calling.bam.bai" for c in "AB"]
    targets.extend([f"trimmed.barcoded.{nr}.fastq.gz" for nr in [1, 2]])
    targets.extend([f"trimmed.non_barcoded.{nr}.fastq.gz" for nr in [1, 2]])
    run(workdir=path, snakemake_args=targets + DEFAULT_SMK_ARGS)
    return path


@pytest.fixture
def workdir(_workdir, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir, path)
    return path


def test_init(tmp_path):
    init(tmp_path / "analysis", TESTDATA_DBS_READ1, "dbs")


def test_config(tmp_path):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_DBS_READ1, "dbs")
    change_config(
        workdir / "blr.yaml",
        [("read_mapper", "bwa"),
         ("hard_filters.snps", "\"'QUAL < 15','lowQual'\"")]  # Nested parameter
    )


def test_default_read_mapper(workdir):
    n_input_fastq_reads = 2 * count_fastq_reads(workdir / "trimmed.barcoded.1.fastq.gz")
    n_input_fastq_reads += 2 * count_fastq_reads(workdir / "trimmed.non_barcoded.1.fastq.gz")
    assert 0 < count_bam_alignments(workdir / "initialmapping.bam") <= n_input_fastq_reads


# Bowtie2, bwa and minimap2 all use the same format so only bowtie2 is tested.
def test_trim_dbs_bowtie2(workdir):
    # bowtie2 is the default so we don't need to rerun anything
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", "bowtie2")]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=trimmed + ["--forcerun", "trim"] + DEFAULT_SMK_ARGS)

    assert 0 < count_fastq_reads(workdir / trimmed[0]) <= count_fastq_reads(TESTDATA_DBS_READ1)
    assert 0 < count_fastq_reads(workdir / trimmed[1]) <= count_fastq_reads(TESTDATA_DBS_READ2)


def test_trim_dbs_ema(workdir):
    # Check that all bins exist and have barcodes in groups.
    for nr in range(5):
        bin_file = workdir / "fastq_bins" / f"ema-bin-00{nr}"
        assert bin_file.exists()
        assert fastq_ema_has_barcodes_grouped(bin_file)

    assert 0 < count_fastq_reads(workdir / "trimmed.barcoded.1.fastq.gz") <= count_fastq_reads(TESTDATA_DBS_READ1)
    assert 0 < count_fastq_reads(workdir / "trimmed.barcoded.2.fastq.gz") <= count_fastq_reads(TESTDATA_DBS_READ2)


@pytest.mark.skipif(shutil.which("lariat") is None, reason="Lariat not installed")
def test_trim_dbs_lariat(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", "lariat")]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=trimmed + ["--forcerun", "trim"] + DEFAULT_SMK_ARGS)
    assert 0 < count_lariat_fastq_reads(workdir / trimmed[0]) <= count_fastq_reads(TESTDATA_DBS_READ1)
    assert fastq_lariat_has_barcodes_grouped(workdir / trimmed[0])


@pytest.fixture(scope="module")
def _workdir_tenx(tmp_path_factory):
    """
    This sets up a workdir for tenx analysis
    """
    path = tmp_path_factory.mktemp(basename="analysis-tenx") / "analysis"
    init(path, TESTDATA_TENX_READ1, "10x")
    change_config(
        path / DEFAULT_CONFIG, [
            ("genome_reference", REFERENCE_GENOME),
            ("barcode_whitelist", TESTDATA_TENX_BARCODES),
            ("chunk_size", "50000"),
            ("phasing_contigs", "null"),
            ("heap_space", "1")
        ]
    )
    return path


@pytest.fixture
def workdir_tenx(_workdir_tenx, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir_tenx, path)
    return path


@pytest.mark.parametrize("read_mapper", ["bowtie2", "ema"])
def test_trim_tenx(workdir_tenx, read_mapper):
    workdir = workdir_tenx
    nr_bins = 5
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", read_mapper),
         ("fastq_bins", str(nr_bins))]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=trimmed + DEFAULT_SMK_ARGS)

    for raw, trimmed in zip((TESTDATA_TENX_READ1, TESTDATA_TENX_READ2), trimmed):
        assert 0 < count_fastq_reads(workdir / trimmed) <= count_fastq_reads(raw)


@pytest.fixture(scope="module")
def _workdir_stlfr(tmp_path_factory):
    """
    This sets up a workdir for tenx analysis
    """
    path = tmp_path_factory.mktemp(basename="analysis-stlfr") / "analysis"
    init(path, TESTDATA_STLFR_READ1, "stlfr")
    change_config(
        path / DEFAULT_CONFIG, [
            ("genome_reference", REFERENCE_GENOME),
            ("chunk_size", "50000"),
            ("phasing_contigs", "null"),
            ("heap_space", "1")
        ]
    )
    return path


@pytest.fixture
def workdir_stlfr(_workdir_stlfr, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir_stlfr, path)
    return path


@pytest.mark.parametrize("read_mapper", ["bowtie2", "ema"])
def test_trim_stlfr(workdir_stlfr, read_mapper):
    workdir = workdir_stlfr
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", read_mapper),
         ("fastq_bins", "5")]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=trimmed + DEFAULT_SMK_ARGS)

    # Check that all bins exist and have barcodes in groups.
    if read_mapper == "ema":
        for nr in range(5):
            bin_file = workdir / "fastq_bins" / f"ema-bin-00{nr}"
            assert bin_file.exists()
            assert fastq_ema_has_barcodes_grouped(bin_file)

    for raw, trimmed in zip((TESTDATA_STLFR_READ1, TESTDATA_STLFR_READ2), trimmed):
        assert 0 < count_fastq_reads(workdir / trimmed) <= count_fastq_reads(raw)


@pytest.mark.skipif(shutil.which("lariat") is None, reason="Lariat not installed")
def test_trim_stlfr_lariat(workdir_stlfr):
    workdir = workdir_stlfr
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", "lariat")]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=trimmed + DEFAULT_SMK_ARGS)
    assert 0 < count_lariat_fastq_reads(workdir / trimmed[0]) <= count_fastq_reads(TESTDATA_STLFR_READ1)
    assert fastq_lariat_has_barcodes_grouped(workdir / trimmed[0])


@pytest.fixture(scope="module")
def _workdir_tellseq(tmp_path_factory):
    """
    This sets up a workdir for tenx analysis
    """
    path = tmp_path_factory.mktemp(basename="analysis-tellseq") / "analysis"
    init(path, TESTDATA_TELLSEQ_READ1, "tellseq")
    change_config(
        path / DEFAULT_CONFIG, [
            ("genome_reference", REFERENCE_GENOME),
            ("tellseq_index", TESTDATA_TELLSEQ_INDEX),
            ("chunk_size", "50000"),
            ("phasing_contigs", "null"),
            ("heap_space", "1")
        ]
    )
    return path


@pytest.fixture
def workdir_tellseq(_workdir_tellseq, tmp_path):
    """Make a fresh copy of the prepared analysis directory"""
    path = tmp_path / "analysis"
    shutil.copytree(_workdir_tellseq, path)
    return path


# Trimming is the same for minimap2, bowtie2, and bwa
def test_trim_tellseq_bowtie2(workdir_tellseq):
    workdir = workdir_tellseq
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", "bowtie2")]
    )
    targets = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=targets + DEFAULT_SMK_ARGS)
    for raw, trimmed in zip((TESTDATA_TELLSEQ_READ1, TESTDATA_TELLSEQ_READ2), targets):
        assert 0 < count_fastq_reads(workdir / trimmed) <= count_fastq_reads(raw)


def test_trim_tellseq_ema(workdir_tellseq):
    workdir = workdir_tellseq
    nr_bins = 5
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", "ema"),
         ("fastq_bins", str(nr_bins))]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    trimmed_nobc = ["trimmed.non_barcoded.1.fastq.gz", "trimmed.non_barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=trimmed + trimmed_nobc + DEFAULT_SMK_ARGS)

    # Check that all bins exist and have barcodes in groups.
    for nr in range(nr_bins):
        bin_file = workdir / "fastq_bins" / f"ema-bin-00{nr}"
        assert bin_file.exists()
        assert fastq_ema_has_barcodes_grouped(bin_file)

    for raw, trim, trim_nobc in zip((TESTDATA_TELLSEQ_READ1, TESTDATA_TELLSEQ_READ2), trimmed, trimmed_nobc):
        count_trimmed = count_fastq_reads(workdir / trim) + count_fastq_reads(workdir / trim_nobc)
        assert 0 < count_trimmed <= count_fastq_reads(raw)


@pytest.mark.skipif(shutil.which("lariat") is None, reason="Lariat not installed")
def test_trim_tellseq_lariat(workdir_tellseq):
    workdir = workdir_tellseq
    change_config(
        workdir / DEFAULT_CONFIG,
        [("read_mapper", "lariat")]
    )
    trimmed = ["trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=trimmed + DEFAULT_SMK_ARGS)
    assert 0 < count_lariat_fastq_reads(workdir / trimmed[0]) <= count_fastq_reads(TESTDATA_TELLSEQ_READ1)
    assert fastq_lariat_has_barcodes_grouped(workdir / trimmed[0])


non_default_mappers = ["bwa", "minimap2", "bowtie2"] + (["lariat"] if shutil.which("lariat") else [])


@pytest.mark.parametrize("read_mapper", non_default_mappers)
def test_nondefault_read_mappers(tmp_path, read_mapper):
    workdir = tmp_path / "analysis"
    init(workdir, TESTDATA_DBS_READ1, "dbs")
    change_config(
        workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME),
         ("read_mapper", read_mapper),
         ("phasing_contigs", "null"),
         ("heap_space", "1"),
         ("fastq_bins", "5")]
    )
    targets = ["initialmapping.bam", "trimmed.barcoded.1.fastq.gz", "trimmed.barcoded.2.fastq.gz"]
    run(workdir=workdir, snakemake_args=targets + DEFAULT_SMK_ARGS)
    if read_mapper == "lariat":
        n_input_fastq_reads = 2 * count_lariat_fastq_reads(workdir / "trimmed.barcoded.1.fastq.gz")
    else:
        n_input_fastq_reads = 2 * count_fastq_reads(workdir / "trimmed.barcoded.1.fastq.gz")
    assert 0 < count_bam_alignments(workdir / "initialmapping.bam") <= n_input_fastq_reads


def test_final_compressed_reads_exist(workdir):
    targets = ["reads.1.final.fastq.gz", "reads.2.final.fastq.gz"]
    run(workdir=workdir, snakemake_args=targets + DEFAULT_SMK_ARGS)
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
    target = "chunks/chrA.calling.bam"
    run(workdir=workdir, snakemake_args=[target] + DEFAULT_SMK_ARGS)
    with pysam.AlignmentFile(workdir / target) as af:
        # Ensure that ApplyBQSR was run on the file by inspecting the @PG lines in the header
        bqsr_header_entries = [entry for entry in af.header["PG"] if entry["ID"] == "GATK ApplyBQSR"]
        assert bqsr_header_entries


variant_callers = ["freebayes", "bcftools", "gatk"] + (["deepvariant"] if shutil.which("deepvariant") else [])


@pytest.mark.parametrize("variant_caller", variant_callers)
def test_call_variants(workdir, variant_caller):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null"), ("variant_caller", variant_caller)]
    )
    target = "chunks/chrA.variants.called.vcf"
    run(workdir=workdir, snakemake_args=[target] + DEFAULT_SMK_ARGS)
    assert workdir.joinpath(target).is_file()


def test_filter_variants(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null"), ("filter_variants", "true")]
    )
    target = "chunks/chrA.variants.called.filtered.vcf"
    run(workdir=workdir, snakemake_args=[target] + DEFAULT_SMK_ARGS)
    assert workdir.joinpath(target).is_file()


def test_haplotag(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", REFERENCE_VARIANTS)]
    )
    target = "chunks/chrA.calling.phased.bam"
    run(workdir=workdir, snakemake_args=[target] + DEFAULT_SMK_ARGS)
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
    run(workdir=old_workdir, snakemake_args=["final.bam", "final.molecule_stats.filtered.tsv"] + DEFAULT_SMK_ARGS)

    # Initialize new dir based on old and run setup.
    init_from_dir(new_workdir, [old_workdir], "dbs")
    change_config(
        new_workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME),
         ("chunk_size", "50000"),
         ("phasing_contigs", "null"),
         ("skip_contigs", "null")]
        )
    run(workdir=new_workdir, snakefile="run_anew.smk", snakemake_args=DEFAULT_SMK_ARGS)


def test_merge_workdirs(tmp_path, workdir):
    other_workdir = tmp_path / "other"
    merge_workdir = tmp_path / "merge"

    # Generate all require files in workdir
    run(workdir=workdir, snakemake_args=["final.bam", "final.molecule_stats.filtered.tsv"] + DEFAULT_SMK_ARGS)

    # Copy and change sample_nr to simulate other library
    shutil.copytree(workdir, other_workdir)
    change_config(
        other_workdir / CONFIGURATION_FILE_NAME,
        [("sample_nr", "2")]
    )

    # Initialize new dir based on old and run setup.
    init_from_dir(merge_workdir, [workdir, other_workdir], "dbs")
    change_config(
        merge_workdir / DEFAULT_CONFIG,
        [("genome_reference", REFERENCE_GENOME),
         ("chunk_size", "50000"),
         ("phasing_contigs", "null"),
         ("skip_contigs", "null")]
        )
    run(workdir=merge_workdir, snakefile="run_anew.smk", snakemake_args=DEFAULT_SMK_ARGS)


def test_lsv_calling(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("reference_variants", "null")]
    )
    target = "chunks/chrA.naibr_sv_calls.tsv"
    run(workdir=workdir, snakemake_args=[target] + DEFAULT_SMK_ARGS)
    assert workdir.joinpath(target).is_file()


def test_phasing_contigs(workdir):
    change_config(
        workdir / DEFAULT_CONFIG,
        [("phasing_contigs", "chrA"),
         ("reference_variants", REFERENCE_VARIANTS)]
    )
    targets = ["final.phased.vcf.gz", "final.phased.vcf.gz.tbi"]
    run(workdir=workdir, snakemake_args=targets + DEFAULT_SMK_ARGS)
    assert chromsome_phased_in_vcf(workdir.joinpath(targets[0]), chromosome="chrA")
    assert not chromsome_phased_in_vcf(workdir.joinpath(targets[0]), chromosome="chrB")
    assert not chromsome_phased_in_vcf(workdir.joinpath(targets[0]), chromosome="chrC")
    assert not chromsome_phased_in_vcf(workdir.joinpath(targets[0]), chromosome="chrD")


def test_multiqc_report_complete(workdir):
    expected_sections = {
        "HapCUT2",
        "Stats",
        "Whatshap",
        "mosdepth",
        "Picard",
        "Cutadapt",
        "FastQC",
        "Samtools"
    }
    targets = ["multiqc_report.html", "multiqc_data"]
    run(workdir=workdir, snakemake_args=targets + DEFAULT_SMK_ARGS)
    assert all((workdir / t).exists() for t in targets)

    sections = set()
    with (workdir / "multiqc_data" / "multiqc_sources.txt").open() as f:
        _ = next(f)  # skip header
        for line in f:
            sections.add(line.split("\t")[0])

    assert sections.issubset(expected_sections)
