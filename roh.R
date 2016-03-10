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
head(runs %>% dplyr::filter(((stops-starts) > 50) & values==1) %>% arrange(-lengths))

#CCGO59$ROH <- ''

