#!/bin/bash
set -xeuo pipefail

uname
samtools --version
bowtie2 --version
minimap2 --version
cutadapt --version
starcode --version
snakemake --version
blr --version
picard MarkDuplicates --version || true
ema

if [[ $(uname) == Darwin ]]; then
  md5="md5 -r"
else
  md5=md5sum
fi

( cd testdata && bwa index chr1mini.fasta )
( cd testdata && bowtie2-build chr1mini.fasta chr1mini.fasta > /dev/null )
( cd testdata && samtools faidx chr1mini.fasta )
if test ! -f testdata/chr1mini.dict; then ( cd testdata && gatk CreateSequenceDictionary -R chr1mini.fasta ); fi

pytest -v tests/

# Test full run on BLR library.
rm -rf outdir-bowtie2
blr init --r1=testdata/blr_reads.1.fastq.gz -l blr outdir-bowtie2
blr config \
    --file outdir-bowtie2/blr.yaml \
    --set genome_reference ../testdata/chr1mini.fasta \
    --set dbSNP ../testdata/dbSNP.chr1mini.vcf.gz \
    --set reference_variants ../testdata/HG002_GRCh38_GIAB_highconf.chr1mini.vcf \
    --set phasing_ground_truth ../testdata/HG002_GRCh38_GIAB_highconf_triophased.chr1mini.vcf \
    --set max_molecules_per_bc 1 \
    --set heap_space 1

pushd outdir-bowtie2
blr run
m=$(samtools view mapped.sorted.tag.bcmerge.mkdup.mol.filt.bam | $md5 | cut -f1 -d" ")
test $m == 61b23ba7b0e7a00f788033729de6bdca

# Cut away columns 2 and 3 as these change order between linux and osx
m2=$(cut -f1,4- mapped.sorted.tag.bcmerge.mkdup.mol.filt.phase | $md5 | cut -f1 -d" ")
test $m2 == e6c512bfeb7cb1b230a9320f22c32937
