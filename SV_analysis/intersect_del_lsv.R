library(InteractionSet)
library(GenomicRanges)
library(bedr)

setwd(getwd())

list.files(getwd(),pattern = "*bedpe")

#files for testing
#M10_merged <- read.table('201026.BLR.M10_merge_final.phased.bam.filtered_large_svcalls.bedpe')
#M10_2_3 <- read.table("201028.BLR.M10.2-3.merge_final.phased.bam.filtered_large_svcalls.bedpe")
#SUSHI_M10_1_2 <- read.table("201119.BLR.SUSHI_M10.1-2.merge_final.phased.bam.filtered_large_svcalls.bedpe")
#M10_2_3_SM10_1_3 <- read.table("201203.BLR.M10.2-3_SM10.1-3_final.phased.bam.filtered_large_svcalls.bedpe")


clean_bedpe <- function(bedpelist) {
  bedpelist <- bedpelist[,-c(3,4,6,12)] # 12 can be kept
  bedpelist <- bedpelist[bedpelist$V7 == "DEL",]
  bedpelist <- makeGRangesFromDataFrame(bedpelist, seqnames.field = 'V1',
                           start.field = 'V2', end.field = 'V5',
                           keep.extra.columns = T) 
  return(bedpelist)
}


M10_merged <- clean_bedpe(M10_merged)
M10_2_3 <- clean_bedpe(M10_2_3)
SUSHI_M10_1_2 <- clean_bedpe(SUSHI_M10_1_2)
M10_2_3_SM10_1_3 <- clean_bedpe(M10_2_3_SM10_1_3)



unique_range1 <- function(range1, range2) {
  remove_lines <- queryHits(findOverlaps(query = range1, subject = range2))
  range1 <- range1[-remove_lines]
  return(range1)
    }

unique_range2 <- function(range2, range1) {
  remove_lines <- queryHits(findOverlaps(query = range1, subject = range2))
  range1 <- range1[-remove_lines]
  return(range1)
}

unique_range1(M10_merged, M10_2_3 )

length(M10_merged)
length(M10_2_3)

length(M10_2_3_SM10_1_3)
length(SUSHI_M10_1_2)

length(unique_range1(M10_merged, M10_2_3 ) )
length(unique_range2(M10_merged, M10_2_3 ) )

length(unique_range1(SUSHI_M10_1_2, M10_2_3_SM10_1_3 ) )

write.table(tumor,file = 'uniqe_tumor_Sdeletions.bed', sep = '\t', quote = F, 
            col.names = F, row.names = F)


# 
# summary(tumor$V10)
# hist(tumor$V10,breaks = 10, xlim=c(20, 130), labels = T)
# 


grange_to_bedpe <- function(df) {
  
  df <- data.frame(df)
  df$end1 <- df$start + 1
  df$end2 <- df$end + 1
  
  df <- df [, c ('seqnames', 'start', 'end1', "seqnames", "end", 'end2',
                 "V7"   ,     "V8"    ,    "V9"  ,     
                 "V10"  ,     "V11"   ,    "V12"  )]
  
  return(df)
}


df <- grange_to_bedpe(tumor)

write.table(df,file = 'Lnew.bedpe', sep = '\t', quote = F, 
            col.names = F, row.names = F)


