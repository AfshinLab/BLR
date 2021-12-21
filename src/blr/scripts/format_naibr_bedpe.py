"""
Fromat BEDPE from NAIBR output TSV

Use indepentant from snakemake with command:

    python format_naibr_bedpe.py input.tsv output.bedpe

"""
import sys

from blr.utils import parse_naibr_tsv


def main(tsv, bedpe):
    with open(tsv) as infile, open(bedpe, "w") as outfile:
        # Header
        print("#chrom1", "start1", "stop1", "chrom2", "start2", "stop2", "sv_type", "sv_id", "sv_length", "qual_score",
              "filter", "info", sep="\t", file=outfile)
        for nr, sv in enumerate(parse_naibr_tsv(infile)):
            # Use BEDPE format according to LinkedSV: https://github.com/WGLab/LinkedSV#-sv-call-file
            print(
                f"chr{sv.chr1}",
                sv.break1,
                sv.break1+1,
                f"chr{sv.chr2}",
                sv.break2,
                sv.break2+1,
                sv.svtype(),
                f"ID{nr:04}",
                len(sv),
                sv.score,
                sv.pass_filter,
                "ZYGOSITY={};NUM_SPLIT_MOLECULES={};NUM_DISCORDANT_READS={};ORIENTATION={};HAPLOTYPE={}".format(
                    sv.zygosity(),
                    sv.nr_split_molecules,
                    sv.nr_discordant_reads,
                    sv.orientation,
                    sv.haplotype_string,
                ),
                sep="\t",
                file=outfile
            )


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if len(sys.argv) == 3:
            tsv = sys.argv[1]
            bedpe = sys.argv[2]
        else:
            print(__doc__)
            sys.exit(1)
    else:
        tsv = snakemake.input.tsv  # noqa: F821
        bedpe = snakemake.output.bedpe  # noqa: F821

    main(tsv, bedpe)
