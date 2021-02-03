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
echo  ${file::${#file}-6}

Rscript filter_naibr.R -f $file -t DEL -q 0.5 #DUP INV  

done


ls *sorted_merged.vcf > files_to_merge.txt

SURVIVOR merge files_to_merge.txt 5000 1 0 0 0 10000 4_naibrs_merged.vcf

Rscript plot_merged_vcf.R -i 4_naibrs_merged.vcf

rm *log
