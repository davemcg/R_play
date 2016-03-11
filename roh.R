library(accelerometry)
library(ggplot2)
library(data.table)
library(dplyr)

# takes output from bcftools roh
CCGO59 <- fread('~/Dropbox/CCGO_800059.bcftools.roh',skip=3)
CCGO60 <- fread('~/Dropbox/CCGO_800060.bcftools.roh',skip=3)
CCGO61 <- fread('~/Dropbox/CCGO_800061.bcftools.roh',skip=3)
CCGO62 <- fread('~/Dropbox/CCGO_800062.bcftools.roh',skip=3)

blockFinder <- function(df, Subject) {
  the_df <- df

  colnames(the_df) <- c("Chr","Pos","HomozygousState","Quality")
  runs <- rle2(the_df$HomozygousState,indices = TRUE)
  runs <- data.table(runs)
  nrow(runs %>% dplyr::filter(((stops-starts) > 50) & values==1) )
  # find roh block with >50 consecutive calls of homozygosity
  over50 <- (runs %>% dplyr::filter(((stops-starts) > 50) & values==1) %>% arrange(-lengths))
  # calculate genomic space in the roh
  over50$GenomicDistance<-apply(over50,1,function(x) abs(the_df[x[2]]$Pos - the_df[x[3]]$Pos))
  # only keep roh with >1,000,000kb
  over50 <- over50 %>% filter(GenomicDistance > 1000000)

  # label original DF with roh blocks
  the_df$ROH <- "No"
  for (i in 1:nrow(over50)){
    the_df[seq(over50[i,]$starts,over50[i,]$stops),"ROH"] <- Subject
  }
  return(the_df)
}
one<-blockFinder(CCGO59,"59")
two<-blockFinder(CCGO60,"60")
three<-blockFinder(CCGO61,"61")
four<-blockFinder(CCGO62,"62")

test <- rbind(one,two,three,four)
test <- test[!duplicated(test),]
test$Chr <- factor(test$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))

ggplot() + geom_point(data=test,aes(x=Pos,y=ROH)) + facet_wrap(~Chr,ncol=2)

var <- cbind(c("chr12","chr7","chr1","chr15"),c(9465464,65548068,153314125,72190416),c(0,0,0,0),c('..','..','..','..'),c(62,62,62,62),c("RP11-22B23.1","ASL","PGLYRP4","MYO9A"))
colnames(var) <- c("Chr","Pos","HomozygousState","Quality","ROH","Gene")
var <- data.table(var)
var$Chr <- factor(var$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
var$Pos <- as.numeric(var$Pos)

test$Gene <- ''
test2<-rbind(test,var)
test3 <- test2 %>% dplyr::filter(ROH != "No")
