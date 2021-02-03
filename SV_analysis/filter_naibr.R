###
# This code is to filter bedpe file generated from Naibr and make it
# compatible with SURVIVOR, with additional filtration of lengths
# and quality. 
###

library(optparse)

system("export LANG=en_US.UTF-8 && export LC_ALL=en_US.UTF-8")

option_list = list(
  make_option(c("-f", "--file"), type="character", default="", 
              help="bedpe file generated from NAIBR", metavar="character"),
  make_option(c("-t", "--type"), type="character", default="DEL", 
              help="Choose which SV type to keep [DEL, DUP, INV] [default= %default]", metavar="character"),
  make_option(c("-M", "--maxlength"), type="integer", default=10000000, 
              help="Max SV length to keep [default= %default]", metavar="integer"),
  make_option(c("-m", "--minlength"), type="integer", default=10000, 
              help="Min SV length to keep [default= %default]", metavar="integer"),
  make_option(c("-q", "--quantile"), type="double", default=0.1, 
              help="Remove SV with QC less than the quantile [default= %default]", metavar="double"),
  make_option(c("-d", "--max_distance"), type="integer", default=1000,
              help="Max distance between SV to be merged [default= %default]", metavar="integer")
); 


##################

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (nchar(opt$file)>0) {
  df <- read.table(opt$file)
  file_name <- opt$file
  file_name <- gsub("\\.[[:alnum:]]+$", "\\1", file_name) 
} else {
  stop("Date parameter must be provided. See script usage (--help)")
}

type <- opt$type
maxlength <- opt$maxlength
minlength <- opt$minlength
quant <- opt$quantile
max_distance <- opt$max_distance

# ## for testing ##
 # df <- read.table('201026.BLR.M10_merge.ema_final.naibr_sv_calls.bedpe')
 # type <- 'DEL'
 # file_name <- 'test'


naibr_to_linkedSV <- function(df) {
  type <- gsub(pattern = "Type=(.*);Zygosity.*", "\\1", x = df$V12)
  col12 <- gsub(pattern = "Type=*(.*?)\\s*;", replacement = "", x = df$V12)
  df <- data.frame(df[,1:6],type,df[,7], (df$V6-df$V2), df[,8],"PASS",col12)
  colnames(df) <- c("#chrom1",	"start1",	"stop1",	"chrom2",	"start2",	"stop2",	"sv_type",
                    "sv_id",	"sv_length",	"qual_score",	"filter",	"info")
  return(df)
}

bedpe_to_bed <- function(df) {
  ## only for deletions, inversion and duplication 
  df <- data.frame(df[,1:2],df[,5],df[,7:ncol(df)])
 return(df)
  
}
bed_to_bedpe <- function(df) {
  ## For plotting in IGV, start1 = end1, start2 = end2, chr1 = chr2
  df <- data.frame(df[,1:2],df[,2],df[,1],df[,3],df[,3],df[,4:ncol(df)])
  return(df)
}


df <- naibr_to_linkedSV(df)
df <- df[df$sv_type==type,]

chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM","chrEBV")
df$`#chrom1` <- factor(df$`#chrom1`, levels = chrOrder) 
df <- df[order(df$`#chrom1`, df$start1),]
df <- bedpe_to_bed(df)
df <- df[df$sv_length <= maxlength,]
df <- df[df$sv_length >= minlength,]

#summary(df)
#hist(df$qual_score,breaks = 100,freq = F)
#length(df$qual_score[df$qual_score<=quant10])
quant <- quantile(df$qual_score,c(quant))
df <- df[df$qual_score>=quant,]

bedfile       <- paste0(file_name,"_filtered_sorted.bed")
bedfile_merged <- paste0(file_name,"_filtered_sorted_merged.bed")
bedpefile_merged   <- paste0(file_name,"_filtered_sorted_merged.bedpe")
vcffile_merged <-  paste0(file_name,"_filtered_sorted_merged.vcf")

write.table(df, file = bedfile ,sep = '\t', row.names = F, col.names = F, quote = F)
system( sprintf("bedtools merge -i %s -c 4,5,6,7,8,9 -o first,collapse,max,mean,first,count -d %d > tmp_%s", bedfile,max_distance,bedfile_merged) )

#removing blacklists and gaps
system( sprintf("bedtools intersect -v -a tmp_%1$s -b hg38_black_list.bed -b hg38_gap.bed > %1$s", bedfile_merged))

dd <- read.table(bedfile_merged)
write.table(bed_to_bedpe(dd), file = bedpefile_merged ,sep = '\t', row.names = F, col.names = F, quote = F)

system( sprintf("SURVIVOR bedpetovcf %s tmp_%s",bedpefile_merged,vcffile_merged))
system( sprintf("less tmp_%1$s | sed s'/STRANDS=[0-9][- 0-9];//'g > %1$s",vcffile_merged))
system( sprintf("rm %s %s tmp_%s tmp_%s",bedfile,bedfile_merged,bedfile_merged,vcffile_merged))
