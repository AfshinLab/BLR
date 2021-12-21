"""
Convert NAIBR output TSV to VCF format

Use indepentant from snakemake with command:

    python format_naibr_vcf.py input.tsv reference.fasta.fai output.vcf
"""
import sys

from blr.utils import parse_fai, parse_naibr_tsv


def main(inputname, fainame, outputname):
    with open(fainame) as f:
        chromosomes = parse_fai(f)

    rename_chroms = any(c.name.startswith("chr") for c in chromosomes)

    # Build header
    header_string = '''##fileformat=VCFv4.2
##source=NAIBR_bedpe_conversion'''

    header_string += ''.join([f"\n##contig=<ID={c.name},length={c.length}>" for c in chromosomes])

    header_string += '''
##FILTER=<ID=PASS,Description="Passed the software filter">
##FILTER=<ID=FAIL,Description="Failed the software filter">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, DUP=Duplication, INV=Inversion">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVSCORE,Number=1,Type=Float,Description="NAIBR Score for SV">
##INFO=<ID=NRREADS,Number=1,Type=Integer,Description="Number of supporting discordant reads">
##INFO=<ID=NRSPLIT,Number=1,Type=Integer,Description="Number of supporting split molecules">
##INFO=<ID=ZYGOSITY,Number=1,Type=String,Description="Zygosity">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE'''

    # Lines list. Needed to output the lines in numerical order at the end
    output = []

    # Reading NAIBR bedpe file and parsing each line
    with open(inputname) as f:

        for nr, sv in enumerate(parse_naibr_tsv(f)):
            # Construction of some fields like ID and INFO
            name = f'NAIBR_{nr}'

            output_infos = ';'.join((
                f'END={sv.break2}',
                f'SVTYPE={sv.svtype()}',
                f'SVLEN={len(sv)}',
                f'SVSCORE={sv.score}',
                f'NRREADS={sv.nr_discordant_reads}',
                f'NRSPLIT={sv.nr_split_molecules}',
                f'ZYGOSITY={sv.zygosity()}',
            ))

            # Final line construction and appending to list
            output_line = (
                f"chr{sv.chr1}" if rename_chroms else sv.chr1,
                sv.break1,
                name,
                'N',
                f'<{sv.svtype()}>',
                sv.score,
                sv.pass_filter,
                output_infos,
                'GT',
                './.'
            )
            output.append(output_line)

    # Final sort and print
    with open(outputname, "w") as f:
        print(header_string, file=f)
        output.sort(key=lambda x: (x[0], x[1]))
        for line in output:
            print(*line, sep="\t", file=f)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if len(sys.argv) == 4:
            tsv = sys.argv[1]
            fai = sys.argv[2]
            vcf = sys.argv[3]
        else:
            print(__doc__)
            sys.exit(1)
    else:
        tsv = snakemake.input.tsv  # noqa: F821
        fai = snakemake.input.fai  # noqa: F821
        vcf = snakemake.output.bedpe  # noqa: F821

    main(tsv, fai, vcf)
