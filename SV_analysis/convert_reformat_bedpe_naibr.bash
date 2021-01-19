#!/bin/bash

sort_vcf () {
local file=$1   
echo 'sorting .. '$file
cat $file | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ${file::${#file}-4}'_sorted.vcf'
gzip_index  ${file::${#file}-4}'_sorted.vcf';
}

gzip_index () {
local file=$1   
echo 'Zipping and indexing .. '$file
bgzip -c $file > ${file}".gz"
tabix ${file}".gz"  
}


################

for file in $(ls *naibr_sv_calls.bedpe)
do
echo  ${file::${#file}-6}"_LSV.bedpe"

Rscript reformat_naibr_output.R -f $file -o ${file::${#file}-6}"_LSV.bedpe"
SURVIVOR bedpetovcf ${file::${#file}-6}"_LSV.bedpe" ${file::${#file}-6}"_LSV.vcf"

less ${file::${#file}-6}"_LSV.vcf" | sed s'/STRANDS=[0-9][- 0-9];//'g >  ${file::${#file}-6}"_LSV_noSTRANDS.vcf"

sort_vcf ${file::${#file}-6}"_LSV_noSTRANDS.vcf"


rm ${file::${#file}-6}"_LSV.vcf"
rm ${file::${#file}-6}"_LSV_noSTRANDS.vcf"
done
