from pathlib import Path
from shutil import copyfile
import subprocess

TESTDATA_STATS = Path("tests/testdata_multiqc_blr/data/example_stats.log")
TESTDATA_STATS_PHASEBLOCK_DATA = Path("tests/testdata_multiqc_blr/data/example_stats_phaseblock_data.tsv")
TESTDATA_STATS_MOLECULES = Path("tests/testdata_multiqc_blr/data/example_stats_molecule_lengths.tsv")
TESTDATA_HAPCUT2_PHASING_STATS = Path("tests/testdata_multiqc_blr/data/example_hapcut2_phasing_stats.txt")
TESTDATA_WHATSHAP_STATS = Path("tests/testdata_multiqc_blr/data/example_whatshap_stats.tsv")


REF_STATS = Path("tests/testdata_multiqc_blr/reference/example_stats.txt")
REF_STATS_PHASEBLOCK_DATA = Path("tests/testdata_multiqc_blr/reference/stats_phaseblock_lengths.txt")
REF_STATS_MOLECULES = Path("tests/testdata_multiqc_blr/reference/stats_molecule_lengths.txt")
REF_HAPCUT2_PHASING_STATS = Path("tests/testdata_multiqc_blr/reference/hapcut2_phasing_stats.txt")
REF_WHATSHAP_STATS = Path("tests/testdata_multiqc_blr/reference/whatshap_stats.txt")
REF_WHATSHAP_STATS_SNVS_PHASED = Path("tests/testdata_multiqc_blr/reference/whatshap_stats_snvs_phased.txt")


def comp_files_linewise(file1: Path, file2: Path):
    with open(file1) as f1, open(file2) as f2:
        for l_f1, l_f2 in zip(f1, f2):
            assert l_f1 == l_f2


def test_stats(tmpdir):
    copyfile(TESTDATA_STATS, tmpdir / "example.log")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    assert Path(tmpdir / "multiqc_data" / "example_stats.txt").exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / "example_stats.txt"),
                        REF_STATS)


def test_stats_phaseblock_data(tmpdir):
    copyfile(TESTDATA_STATS_PHASEBLOCK_DATA, tmpdir / "example.phaseblock_data.tsv")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    assert Path(tmpdir / "multiqc_data" / "stats_phaseblock_lengths.txt").exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / "stats_phaseblock_lengths.txt"),
                        REF_STATS_PHASEBLOCK_DATA)


def test_stats_molecule_lengths(tmpdir):
    copyfile(TESTDATA_STATS_MOLECULES, tmpdir / "example.molecule_lengths.tsv")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "stats"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    assert Path(tmpdir / "multiqc_data" / "stats_molecule_lengths.txt").exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / "stats_molecule_lengths.txt"),
                        REF_STATS_MOLECULES)


def test_hapcut2(tmpdir):
    copyfile(TESTDATA_HAPCUT2_PHASING_STATS, tmpdir / "example.txt")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "hapcut2"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    assert Path(tmpdir / "multiqc_data" / "hapcut2_phasing_stats.txt").exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / "hapcut2_phasing_stats.txt"),
                        REF_HAPCUT2_PHASING_STATS)


def test_whatshap(tmpdir):
    copyfile(TESTDATA_WHATSHAP_STATS, tmpdir / "example.whatshap_stats.tsv")

    subprocess.run(["multiqc", "-f", tmpdir, "-o", tmpdir, "-m", "whatshap"])

    assert Path(tmpdir / "multiqc_report.html").exists()
    assert Path(tmpdir / "multiqc_data" / "whatshap_stats.txt").exists()
    assert Path(tmpdir / "multiqc_data" / "whatshap_stats_snvs_phased.txt").exists()

    comp_files_linewise(Path(tmpdir / "multiqc_data" / "whatshap_stats.txt"),
                        REF_WHATSHAP_STATS)
    comp_files_linewise(Path(tmpdir / "multiqc_data" / "whatshap_stats_snvs_phased.txt"),
                        REF_WHATSHAP_STATS_SNVS_PHASED)
