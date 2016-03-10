library(accelerometry)
library(ggplot2)
library(data.table)
library(dplyr)

# takes output from bcftools roh
CCGO59 <- fread('~/Dropbox/CCGO_800059.bcftools.roh',skip=3)
CCGO60 <- fread('~/Dropbox/CCGO_800060.bcftools.roh',skip=3)
CCGO61 <- fread('~/Dropbox/CCGO_800061.bcftools.roh',skip=3)
CCGO62 <- fread('~/Dropbox/CCGO_800062.bcftools.roh',skip=3)

the_df <- CCGO62

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

# http://stackoverflow.com/questions/8091303/simultaneously-merge-multiple-data-frames-in-a-list
my.list <- list(one,two,three,four)
# make unique names
my.list2 = Map(function(x, i) setNames(x, ifelse(names(x) %in% match.by,names(x), sprintf('%s.%d', names(x), i))), my.list, seq_along(my.list))
merged.data.frame = Reduce(function(...) merge(..., by=c("Chr","Pos"),all=T,allow.cartesian=TRUE), my.list2)

merged.data.frame$Chr <- sapply(merged.data.frame$Chr,function(x) substr(x,4,6))
ggplot(merged.data.frame) + geom_point(aes(x=Pos,y=ROH.1)) + geom_point(aes(x=Pos,y=ROH.2)) + geom_point(aes(x=Pos,y=ROH.3)) + geom_point(aes(x=Pos,y=ROH.4)) + facet_wrap(~Chr,scales='free_x',ncol=2)