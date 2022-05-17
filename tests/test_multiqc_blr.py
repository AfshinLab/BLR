from pathlib import Path
from shutil import copyfile
import subprocess

TESTDATA_BASE = Path("tests/testdata_multiqc_blr/data")
TESTDATA_STATS = TESTDATA_BASE / "example_stats.log"
TESTDATA_STATS_PHASEBLOCK_DATA = TESTDATA_BASE / "example_stats_phaseblock_data.tsv"
TESTDATA_STATS_MOLECULES = TESTDATA_BASE / "example_stats_molecule_lengths.tsv"
TESTDATA_STATS_MOLECULE_STATS = TESTDATA_BASE / "example.molecule_stats.txt"
TESTDATA_STATS_BARCODE_STATS = TESTDATA_BASE / "example.barcode_stats.txt"
TESTDATA_HAPCUT2_PHASING_STATS = TESTDATA_BASE / "example_hapcut2_phasing_stats.txt"
TESTDATA_HAPCUT2_PHASING_STATS_AUN = TESTDATA_BASE / "example_hapcut2_phasing_stats_auN.txt"
TESTDATA_HAPCUT2_PHASING_STATS_CHROM = TESTDATA_BASE / "example_hapcut2_phasing_stats_chroms.txt"
TESTDATA_WHATSHAP_STATS = TESTDATA_BASE / "example_whatshap_stats.tsv"
TESTDATA_WHATSHAP_HAPLOTAG = TESTDATA_BASE / "example.haplotag.log"

REF_BASE = Path("tests/testdata_multiqc_blr/reference")


def comp_files_linewise(file1: Path, file2: Path):
    with open(file1) as f1, open(file2) as f2:
        for l_f1, l_f2 in zip(f1, f2):
            assert l_f1 == l_f2


def test_stats(tmpdir):
    copyfile(TESTDATA_STATS, tmpdir / "example.log")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    file = "example_stats.txt"
    assert Path(tmpdir / "multiqc_data" / file).exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)


def test_stats_phaseblock_data(tmpdir):
    copyfile(TESTDATA_STATS_PHASEBLOCK_DATA, tmpdir / "example.phaseblock_data.tsv")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    file = "stats_phaseblock_lengths.txt"
    assert Path(tmpdir / "multiqc_data" / file).exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)


def test_stats_molecule_lengths(tmpdir):
    copyfile(TESTDATA_STATS_MOLECULES, tmpdir / "example.molecule_lengths.tsv")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    file = "stats_molecule_lengths.txt"
    assert Path(tmpdir / "multiqc_data" / file).exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)


def test_stats_molecule_stats(tmpdir):
    copyfile(TESTDATA_STATS_MOLECULE_STATS, tmpdir / "example.molecule_stats.txt")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    output_files = [
        "multiqc_molecule_stats_RB.txt",
        "multiqc_molecule_stats_MB.txt",
        "multiqc_molecule_stats_MC.txt",
        "multiqc_molecule_stats_summary.txt",
    ]
    for file in output_files:
        assert Path(tmpdir / "multiqc_data" / file).exists()

        comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)


def test_stats_barcode_stats(tmpdir):
    copyfile(TESTDATA_STATS_BARCODE_STATS, tmpdir / "example.barcode_stats.txt")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    output_files = [
        "multiqc_barcode_stats_RB.txt",
        "multiqc_barcode_stats_CB.txt",
        "multiqc_barcode_stats_summary.txt",
    ]
    for file in output_files:
        assert Path(tmpdir / "multiqc_data" / file).exists()

        comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)


def test_hapcut2(tmpdir):
    copyfile(TESTDATA_HAPCUT2_PHASING_STATS, tmpdir / "example.txt")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "hapcut2"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    file = "hapcut2_phasing_stats.txt"
    assert Path(tmpdir / "multiqc_data" / file).exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)


def test_hapcut2_with_auN(tmpdir):
    copyfile(TESTDATA_HAPCUT2_PHASING_STATS_AUN, tmpdir / "example.txt")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "hapcut2"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    file = "hapcut2_phasing_stats.txt"
    assert Path(tmpdir / "multiqc_data" / file).exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / "hapcut2_phasing_stats_auN.txt")


def test_hapcut2_chroms(tmpdir):
    copyfile(TESTDATA_HAPCUT2_PHASING_STATS_CHROM, tmpdir / "example.txt")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "hapcut2"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    file = "hapcut2_phasing_stats.txt"
    assert Path(tmpdir / "multiqc_data" / file).exists()


def test_whatshap_stats(tmpdir):
    copyfile(TESTDATA_WHATSHAP_STATS, tmpdir / "example.whatshap_stats.tsv")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "whatshap"])

    assert Path(tmpdir / "multiqc_report.html").exists()

    for file in ["whatshap_stats.txt", "whatshap_stats_snvs_phased.txt"]:
        assert Path(tmpdir / "multiqc_data" / file).exists()
        comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)


def test_whatshap_haplotag(tmpdir):
    copyfile(TESTDATA_WHATSHAP_HAPLOTAG, tmpdir / "example.haplotag.log")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "whatshap"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    file = "whatshap_haplotag.txt"
    assert Path(tmpdir / "multiqc_data" / file).exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / file), REF_BASE / file)
