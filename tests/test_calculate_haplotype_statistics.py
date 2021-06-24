
from blr.cli.calculate_haplotype_statistics import error_rate_calc


def test_error_rate_calc():
    query_blocklist = [
        [(1, 10, 0, 1, 'A', 'T', None),
         (2, 20, 0, 1, 'G', 'T', None)],
        [(4, 30, 0, 1, 'G', 'T', None),
         (5, 40, 0, 1, 'C', 'T', None),
         (7, 50, 0, 1, 'G', 'C', None),
         (8, 60, 0, 1, 'C', 'T', None)],
        [(10, 70, 0, 1, 'C', 'A', None),
         (11, 80, 0, 1, 'G', 'A', None),
         (13, 110, 0, 1, 'T', 'A', None)],
        [(15, 120, 0, 1, 'C', 'T', None),
         (16, 130, 0, 1, 'G', 'C', None),
         (18, 150, 0, 1, 'T', 'C', None),
         (20, 190, 0, 1, 'A', 'T', None)],
    ]
    reference_blocklist = [
        [(1, 10, 0, 1, 'A', 'T', None),
         (2, 20, 0, 1, 'G', 'T', None)],
        [(4, 30, 1, 0, 'G', 'T', None),
         (5, 40, 1, 0, 'C', 'T', None),
         (7, 50, 0, 1, 'G', 'C', None),
         (8, 60, 0, 1, 'C', 'T', None)],
        [(10, 70, 0, 1, 'C', 'A', None),
         (11, 80, 1, 0, 'G', 'A', None),
         (13, 110, 0, 1, 'T', 'A', None)],
        [(15, 120, 0, 1, 'C', 'T', None),
         (16, 130, 0, 1, 'G', 'C', None),
         (18, 150, 0, 1, 'T', 'C', None),
         (20, 190, 0, 1, 'A', 'T', None)]
    ]
    error_result = error_rate_calc(query_blocklist, reference_blocklist, "chrA", num_snps=13)
    assert error_result.get_switch_rate() == 2 / 4
    assert error_result.get_mismatch_rate() == 1 / 13
    assert error_result.get_N50_phased_portion() == 40
    assert error_result.get_AN50() == 30.0
