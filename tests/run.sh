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
mosdepth --version

if [[ $(uname) == Darwin ]]; then
  md5="md5 -r"
else
  md5=md5sum
fi

pushd blr-testdata
bwa index ref.fasta
bowtie2-build ref.fasta ref.fasta > /dev/null
samtools faidx ref.fasta
test -f ref.dict || gatk CreateSequenceDictionary -R ref.fasta
popd

pytest -v tests/

# Test full run on BLR library.
rm -rf outdir-bowtie2
blr init --r1=blr-testdata/dbs_reads.1.fastq.gz -l dbs outdir-bowtie2
blr config \
    --file outdir-bowtie2/blr.yaml \
    --set genome_reference ../blr-testdata/ref.fasta \
    --set dbSNP ../blr-testdata/dbSNP.vcf.gz \
    --set reference_variants ../blr-testdata/HG002_GRCh38_GIAB_highconf.vcf.gz \
    --set phasing_ground_truth ../blr-testdata/HG002_GRCh38_GIAB_highconf_triophased.vcf.gz \
    --set max_molecules_per_bc 1 \
    --set heap_space 1 \
    --set chunk_size 10000 \
    --set phasing_contigs null \
    --set contigs_skipped null

cd outdir-bowtie2
blr run
m=$(samtools view final.bam | $md5 | cut -f1 -d" ")
test $m == 73d693cb6ce4f3e17c469670ab5f54a0
