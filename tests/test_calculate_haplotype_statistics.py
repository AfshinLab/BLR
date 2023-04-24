
from blr.cli.calculate_haplotype_statistics import error_rate_calc, get_chrom_lengths_from_vcf


def test_error_rate_calc():
    query_blocklist = [
        [(1, 10, (0, 1), ('A', 'T', None)),
         (2, 20, (0, 1), ('G', 'T', None))],
        [(4, 30, (0, 1), ('G', 'T', None)),
         (5, 40, (0, 1), ('C', 'T', None)),
         (7, 50, (0, 1), ('G', 'C', None)),
         (8, 60, (0, 1), ('C', 'T', None))],
        [(10, 70, (0, 1), ('C', 'A', None)),
         (11, 80, (0, 1), ('G', 'A', None)),
         (13, 110, (0, 1), ('T', 'A', None))],
        [(15, 120, (0, 1), ('C', 'T', None)),
         (16, 130, (0, 1), ('G', 'C', None)),
         (18, 150, (0, 1), ('T', 'C', None)),
         (19, 160, (0, 1), ('A', 'C', None)),
         (20, 190, (0, 1), ('A', 'T', None))],
    ]
    reference_blocklist = [
        [(1, 10, (0, 1), ('A', 'T', None)),
         (2, 20, (0, 1), ('G', 'T', None))],
        [(4, 30, (1, 0), ('G', 'T', None)),
         (5, 40, (1, 0), ('C', 'T', None)),
         (7, 50, (0, 1), ('G', 'C', None)),
         (8, 60, (0, 1), ('C', 'T', None))],
        [(10, 70, (0, 1), ('C', 'A', None)),
         (11, 80, (1, 0), ('G', 'A', None)),
         (13, 110, (0, 1), ('T', 'A', None))],
        [(15, 120, (1, 0), ('C', 'T', None)),
         (16, 130, (1, 0), ('G', 'C', None)),
         (18, 150, (0, 1), ('T', 'C', None)),
         (19, 160, (0, 1), ('A', 'C', None)),
         (20, 190, (1, 0), ('A', 'T', None))]
    ]
    error_result = error_rate_calc(query_blocklist, reference_blocklist, "chrA", num_snps=14)
    assert error_result.get_switch_count() == 2
    assert error_result.get_switch_positions() == 3
    assert error_result.get_switch_rate() == 2 / 3
    assert error_result.get_mismatch_count() == 2
    assert error_result.get_mismatch_positions() == 14
    assert error_result.get_mismatch_rate() == 2 / 14
    assert error_result.get_flat_count() == 5
    assert error_result.get_flat_positions() == 14
    assert error_result.get_flat_error_rate() == 5 / 14
    assert error_result.get_N50() == 40
    assert error_result.get_AN50() == 30.0


def test_get_chrom_lengths_from_vcf(tmp_path):
    vcf = tmp_path / "my.vcf"
    s = "##fileformat=VCFv4.1\n"
    s += "##contig=<ID=A, length=50000>\n"
    s += "##contig=<ID=B, length=30000>\n"
    s += "##contig=<ID=C, length=10000>\n"
    s += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMYSAMPLE\n"
    with open(vcf, "w") as f:
        f.write(s)

    chrom_lengths = get_chrom_lengths_from_vcf(vcf, ["A", "B"])
    assert "A" in chrom_lengths
    assert "B" in chrom_lengths
    assert "C" not in chrom_lengths
    assert chrom_lengths["A"] == 50000
    assert chrom_lengths["B"] == 30000
    assert len(chrom_lengths) == 2
