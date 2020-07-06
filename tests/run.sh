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

pushd blr-testdata
bwa index ref.fasta
bowtie2-build ref.fasta ref.fasta > /dev/null
samtools faidx ref.fasta
test -f ref.dict || gatk CreateSequenceDictionary -R ref.fasta
popd

pytest -v tests/

# Test full run on BLR library.
rm -rf outdir-bowtie2
blr init --r1=blr-testdata/blr_reads.1.fastq.gz -l blr outdir-bowtie2
blr config \
    --file outdir-bowtie2/blr.yaml \
    --set genome_reference ../blr-testdata/ref.fasta \
    --set dbSNP ../blr-testdata/dbSNP.vcf.gz \
    --set reference_variants ../blr-testdata/HG002_GRCh38_GIAB_highconf.vcf \
    --set phasing_ground_truth ../blr-testdata/HG002_GRCh38_GIAB_highconf_triophased.vcf \
    --set max_molecules_per_bc 1 \
    --set heap_space 1 \
    --set chunk_size 10000

cd outdir-bowtie2
blr run
m=$(samtools view final.bam | $md5 | cut -f1 -d" ")
test $m == 26a5649edcc722c9fd762829282fe767

# Cut away columns 2 and 3 as these change order between linux and osx
m2=$(grep -v "^##" final.phased.vcf | $md5 | cut -f1 -d" ")
test $m2 == 4948c00636f686bde2f646c3c7efca72
