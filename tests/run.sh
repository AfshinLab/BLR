#!/bin/bash
set -xeuo pipefail

uname
python --version
samtools --version
bowtie2 --version
minimap2 --version
cutadapt --version
starcode --version
snakemake --version
blr --version
picard MarkDuplicates --version || true
ema
bwa 2>&1 | head -5 || true
mosdepth --version
freebayes --version
vcffilter 2>&1 | head -1 || true
fastqc --version
multiqc --version
bcftools --version
gatk --version
whatshap --version

pushd blr-testdata
bwa index ref.fasta
bowtie2-build ref.fasta ref.fasta > /dev/null
samtools faidx ref.fasta
test -f ref.dict || gatk CreateSequenceDictionary -R ref.fasta
popd

pytest -v tests/
