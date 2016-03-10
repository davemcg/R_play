library(accelerometry)
library(ggplot2)
library(data.table)
library(dplyr)

# takes output from bcftools roh
CCGO59 <- fread('~/Dropbox/CCGO_800059.bcftools.roh',skip=3)
CCGO60 <- fread('~/Dropbox/CCGO_800060.bcftools.roh',skip=3)
CCGO61 <- fread('~/Dropbox/CCGO_800061.bcftools.roh',skip=3)
CCGO62 <- fread('~/Dropbox/CCGO_800062.bcftools.roh',skip=3)

colnames(CCGO59) <- c("Chr","Pos","HomozygousState","Quality")
ggplot(CCGO59,aes(x=Pos,y=HomozygousState)) + geom_point() + facet_wrap(~Chr,scales='free_x',ncol=2)
runs <- rle2(CCGO59$HomozygousState,indices = TRUE)
runs <- data.table(runs)
nrow(runs %>% dplyr::filter(((stops-starts) > 50) & values==1) )
# find roh block with >50 consecutive calls of homozygosity
over50 <- (runs %>% dplyr::filter(((stops-starts) > 50) & values==1) %>% arrange(-lengths))
# calculate genomic space in the roh
over50$GenomicDistance<-apply(over50,1,function(x) abs(CCGO59[x[2]]$Pos - CCGO59[x[3]]$Pos))
# only keep roh with >1,000,000kb
over50 <- over50 %>% filter(GenomicDistance > 1000000)

# label original DF with roh blocks
CCGO59$ROH <- "No"
for (i in 1:nrow(over50)){
  CCGO59[seq(over50[i,]$starts,over50[i,]$stops),"ROH"] <- "Yes"
}

ggplot(CCGO59,aes(x=Pos,y=HomozygousState,colour=ROH)) + geom_point() + facet_wrap(~Chr,scales='free_x',ncol=2)
