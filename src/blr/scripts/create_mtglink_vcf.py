"""
Create a VCF file from the FASTA file outputted by MTG-Link and the Manta VCF file used as input to MTG-Link.
Update variants with the assembled insert from the FASTA.

- If there are multiple solutions for one gap the first one is used.
- if there are muliple variants in one gap, the gap is skipped.

"""
import argparse
from pathlib import Path
import dataclasses
from collections import defaultdict

import pysam


# TODO - This currently works in a very basic manner and needs to be improved. 
#        Ideally we should do a realignment of the assembled sequence to the reference
#        to get a more accurate representation of the insert. 
#        Relevant discussion: https://github.com/anne-gcd/MTG-Link/issues/25 

@dataclasses.dataclass
class Gap:
    name: str
    chrom: str
    start: int
    end: int
    direction: str
    kmer: str
    contig_length: int
    qual: str
    record: pysam.FastxRecord = None


COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def revcomp(seq):
    return "".join(COMPLEMENT.get(base, base) for base in reversed(seq))


def parse_record(record):
    # >chr1_0-13239699-L+:chr1_13239799-27878410-R+_fwd_1/1.k61 _ len.1085_qual.AB 
    record_id = record.name

    info = record.comment.split(" ")[1]
    length = int(info.split("len.")[1].split("_")[0])
    qual = info.split("qual.")[1]

    *gap_name, direction, kmer = record_id.split("_")
    kmer = kmer.split(".")[1]
    gap_name = "_".join(gap_name)
    gap_chr =  gap_name.split(":")[0].split("_")[0]
    gap_start = int(gap_name.split(":")[0].split("-")[1])
    gap_end = int(gap_name.split(":")[1].split("-")[0].split("_")[1])
    return Gap(gap_name, gap_chr, gap_start, gap_end, direction, kmer, length, qual, record)    


def main(fasta, vcf_in, vcf_out, flank, debug=False):
    records = defaultdict(list)
    print("Reading FASTA", fasta)
    with pysam.FastxFile(fasta) as f:
        for record in f:
            gap = parse_record(record)
            records[gap.name].append(gap)

    n_read = 0
    n_written = 0
    print("Reading VCF", vcf_in)
    with pysam.VariantFile(vcf_in) as reader:
        sample = reader.header.samples[0]
        print("Writing VCF", vcf_out)
        with pysam.VariantFile(vcf_out, "w", header=reader.header) as writer:
            for nr, (name, gap) in enumerate(records.items(), start=1):
                n_read += 1
                gap = gap[0]

                contig = gap.record.sequence if gap.direction == "fwd" else revcomp(gap.record.sequence)
                insert = contig[flank:-flank]

                if debug:
                    print(gap.name, gap.chrom, gap.start, gap.end)
                
                variants = list(reader.fetch(gap.chrom, gap.start, gap.end))
                variants = [v for v in variants if v.info["SVTYPE"] == "INS"]
    
                if debug:
                    for nr, v in enumerate(variants):
                        print("VARIANT", nr)
                        manta_left = variant.info["LEFT_SVINSSEQ"][0]
                        manta_left_cap = min(100, len(manta_left), len(insert))
                        manta_right = variant.info["RIGHT_SVINSSEQ"][0]
                        manta_right_cap = min(100, len(manta_right), len(insert))
                        print("MANTA", manta_left[:manta_left_cap])
                        print("MTG  ", insert[:manta_left_cap])
                        print("MANTA", manta_right[-manta_right_cap:])
                        print("MTG  ", insert[-manta_right_cap:])
                    
                    print("-" * 80)

                if len(variants) > 1:
                    print("More than one variant in gap", gap.name)
                    for v in variants:
                        print("-", v)
                    continue
    
                if len(variants) == 0:
                    print("No variant in gap", gap.name)
                    continue
    
                variant = variants[0]
                
                # Update variant with assembled insert sequence
                variant.alts = [variant.ref + insert]
                variant.info["SVLEN"] = len(insert)
                variant.id = f"MTGLINK:{nr}"

                writer.write(variant)
                n_written += 1

    print(f"Read {n_read} gaps, wrote {n_written} variants")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        fasta = snakemake.input.fasta  # noqa: F821
        vcf_in = snakemake.input.vcf  # noqa: F821
        vcf_out = snakemake.output.vcf  # noqa: F821
        flank = snakemake.params.flank  # noqa: F821
        log = snakemake.log[0]  # noqa: F821
        debug = False
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("fasta", type=Path, help="Main FASTA file outputted by MTG-Link e.g '*assembled_sequences.fasta'")
        parser.add_argument("vcf", type=Path, help="VCF used for input to MTG-Link")
        parser.add_argument("-o", "--output", type=Path, help="Output VCF", default="-")
        parser.add_argument(
            "--flank", default=550, type=int, 
            help="Length of sequence flanking the gap of each side. Will be removed from"
                "the conting. Default: %(default)s"
        )
        parser.add_argument("--debug", action="store_true", help="Debug mode")
        args = parser.parse_args()

        fasta = args.fasta
        vcf_in = args.vcf
        vcf_out = args.output
        flank = args.flank
        debug = args.debug
        log = f"{vcf_out}.log"

    # Write stdout to log file
    with open(log, "w") as sys.stdout:
        main(fasta, vcf_in, vcf_out, flank, debug)
