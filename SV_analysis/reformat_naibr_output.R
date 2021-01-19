library("optparse")

system("defaults write org.R-project.R force.LANG en_US.UTF-8")

option_list = list(
  make_option(c("-f", "--file"), type="character", default="file", 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--type"), type="character", default="DEL", 
              help="Choose which SV type to use [DEL, DUP, INV] [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 


##################
# 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);
# 
# df <- read.table(opt$file)
# type <- opt$type


## for testing ##
df <- read.table('201026.BLR.M10_merge.ema_final.naibr_sv_calls.bedpe')
type <- 'DEL'


naibr_to_linkedSV <- function(df) {
  
  type <- gsub(pattern = "Type=(.*);Zygosity.*", "\\1", x = df$V12)
  
  col12 <- gsub(pattern = "Type=*(.*?)\\s*;", replacement = "", x = df$V12)
  
  df <- data.frame(df[,1:6],type,df[,7], (df$V6-df$V2), df[,8],"PASS",col12)
  
  colnames(df) <- c("#chrom1",	"start1",	"stop1",	"chrom2",	"start2",	"stop2",	"sv_type",
                    "sv_id",	"sv_length",	"qual_score",	"filter",	"info")
  
  return(df)
}

df <- naibr_to_linkedSV(df)
#df <- df[df$sv_type=="DEL",]

## We need to merge SVs if:
##  - they are of the same type (deletions duplications, and inversions)
## AND
##  -- start are within xx bp
##  -- end within xx bp
##  OR
##  -- they overlap in more than 90 
# 
# library(HelloRanges)
# 
# testpair <- Pairs(first = df,second = df,names =)
# findOverlapPairs(query = df,subject = df)

merge_two_lines <- function(line1,line2) {
  if (line1$sv_type==line2$sv_type && line1[1]==line2[1] && line1[4]==line2[4]) {
    line3 <- line1
    line3$start1 <- min(line1$start1, line2$start1)
    line3$stop1  <- min(line1$stop1 , line2$stop1)
    
    line3$start2 <- max(line1$start2, line2$start2)
    line3$stop2  <- max(line1$stop2 , line2$stop2)
    
  line3$sv_length <- line3$stop2  - line3$start1
  line3$sv_id <- paste0(line1$sv_id,'_',line2$sv_id)
  line3$info  <- paste0(line1$info,'_',line2$info)
  return(line3)
  } else {return(rbind(line1,line2))}
}



library(GenomicRanges)
check_intersect <- function(line1,line2) {
 range1 <- GRanges ( paste0(line1[1],":",line1$start1,"-",line1$stop2) )
 range2 <- GRanges ( paste0(line2[1],":",line2$start1,"-",line2$stop2) )
 bins1 <- disjoin(range1)
 bins2 <- disjoin(range2)

 number <- countOverlaps(bins1, bins2)
 return(number)
      #VERY SLOW
      # x1 <- as.numeric(line1[2])
      # x2 <- as.numeric(line1[5])
      # x3 <- as.numeric(line2[2])
      # x4 <- as.numeric(line2[5])
      # a <- length( intersect(seq(x1,x2,1), seq(x3,x4,1)) )
      # return(a)
  }


# check_intersect(df[6,],df[4,])
# 
# merge_two_lines(df[6,],df[7,])
# 
# df <- df[1:1000,]
# 
# df2 <- data.frame()
# for (row1 in 1:nrow(df)) {
#   for (row2 in row1:nrow(df)) {
#     line1 <- df[row1,]
#     line2 <- df[row2,]
#     if (line1==line2) next
#     if(check_intersect(line1,line2))
#        {df2 <- rbind(df2,merge_two_lines(line1,line2))}
#   }
# }


write.table(file = opt$out, naibr_to_linkedSV(df),
            quote = F,row.names = F,col.names = F,sep = "\t")
