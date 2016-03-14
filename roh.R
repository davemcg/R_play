library(accelerometry)
library(ggplot2)
library(data.table)
library(dplyr)

# takes output from bcftools roh
CCGO59 <- fread('~/Dropbox/CCGO_800059.bcftools.roh2',skip=3)
CCGO60 <- fread('~/Dropbox/CCGO_800060.bcftools.roh2',skip=3)
CCGO61 <- fread('~/Dropbox/CCGO_800061.bcftools.roh2',skip=3)
CCGO62 <- fread('~/Dropbox/CCGO_800062.bcftools.roh2',skip=3)

CCGO159 <- fread('~/Dropbox/CCGO_800159.bcftools.roh',skip=3)
CCGO160 <- fread('~/Dropbox/CCGO_800160.bcftools.roh',skip=3)
CCGO161 <- fread('~/Dropbox/CCGO_800161.bcftools.roh',skip=3)


blockFinder <- function(df, Subject) {
  the_df <- df

  colnames(the_df) <- c("Chr","Pos","HomozygousState","Quality")
  runs <- rle2(the_df$HomozygousState,indices = TRUE)
  runs <- data.table(runs)
  # find roh block with >n (100) consecutive calls of homozygosity
  over50 <- (runs %>% dplyr::filter(((stops-starts) > 100) & values==1) %>% arrange(-lengths))
  # calculate genomic space in the roh
  over50$GenomicDistance <- apply(over50,1,function(x) abs(the_df[x[2]]$Pos - the_df[x[3]]$Pos))
  # only keep roh with >1,000,000kb
  over50 <- over50 %>% filter(GenomicDistance > 1000000)

  # label original DF with roh blocks
  the_df$ROH <- "Coverage"
  for (i in 1:nrow(over50)){
    the_df[seq(over50[i,]$starts,over50[i,]$stops),"ROH"] <- Subject
  }
  return(the_df)
}

blockFinderAppend <- function(df, Subject) {
  the_df <- df
  colnames(the_df) <- c("Chr","Pos","HomozygousState","Quality")
  the_df$End <- the_df$Pos 
  the_df$Subject <- Subject
  the_df$Class <- "Variant"
  the_df <- the_df %>% dplyr::select(Chr, Pos, End, HomozygousState, Quality, Subject, Class)

  # do rle by chr. Still not working????
  chroms <- unique(the_df$Chr)
  # setup block df for ROH blocks
  blocks <- data.table(rbind(c(0,0,0,0,0,0,0)))
  colnames(blocks) <- c('Chr', 'Pos', 'End', 'HomozygousState', 'Quality', 'Subject', 'Class')
  for (i in chroms){
    chr_df <- subset(the_df,Chr==i)
    chr_runs <- data.table(rle2(chr_df$HomozygousState,indices=TRUE))
    over <- chr_runs %>% dplyr::filter(((stops-starts) > 100) & values==1) %>% arrange(-lengths)
    # calculate genomic space in the roh
    over$GenomicDistance<-apply(over,1,function(x) abs(chr_df[x[2]]$Pos - chr_df[x[3]]$Pos))
    over$Pos <- apply(over,1,function(x) chr_df[x[2]]$Pos-0)
    over$End <- apply(over,1,function(x) chr_df[x[3]]$Pos-0)
    over$Chr <- i
    # only keep roh with >1,000,000kb
    over <- over %>% filter(GenomicDistance > 1000000)
    over$HomozygousState <- "1"
    over$Quality <- "0"
    over$Subject <- Subject
    over$Class <- "Block"
    over <- over %>% select(Chr, Pos, End, HomozygousState, Quality, Subject, Class)
    blocks <- rbind(blocks,over)
  }
  output <- rbind(the_df, blocks)
  output$Pos <- as.numeric(output$Pos)
  output$End <- as.numeric(output$End)
  output$HomozygousState <- as.numeric(output$HomozygousState)
  return(output)
}


test59.62 <- rbind(blockFinder(CCGO59,"59"),blockFinder(CCGO60,"60"),blockFinder(CCGO61,"61"),blockFinder(CCGO62,"62"))
# and/or
test159.161 <- rbind(blockFinder(CCGO159,"159"),blockFinder(CCGO160,"160"),blockFinder(CCGO161,"161"))

test59.62$Family<-"CCGO_800059-62"
test159.161$Family<-"CCGO_800159-161"
test <- rbind(test59.62,test159.161)


test <- test[!duplicated(test),]
test$Chr <- factor(test$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))

ggplot() + geom_point(data=test,aes(x=Pos,y=ROH,colour=Family),size=0.5) + facet_wrap(~Chr,ncol=2) + theme_bw() + xlab("") + ylab("Subject")



#### play stuff

# hom rec variants for CCGO_800059-62
var <- cbind(c("chr12","chr7","chr1","chr15"),c(9465464,65548068,153314125,72190416),c(0,0,0,0),c('..','..','..','..'),c(62,62,62,62),rep("CCGO_800059-62.var"),c("RP11-22B23.1","ASL","PGLYRP4","MYO9A"))
colnames(var) <- c("Chr","Pos","HomozygousState","Quality","ROH","Family","Gene")
var <- data.table(var)
var$Chr <- factor(var$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
var$Pos <- as.numeric(var$Pos)
var$ROH <- "62.var"

test$Gene <- ''
test2<-rbind(test,var)
test3 <- test2 %>% dplyr::filter(ROH != "Coverage")
ggplot(data=test3,aes(x=Pos,y=ROH,colour=Family,shape=Gene)) + geom_point(size=0.3) + facet_wrap(~Chr,ncol=3) + theme_bw()

four$Chr <- factor(four$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
ggplot(data=one,aes(x=Pos,y=HomozygousState,colour=ROH,shape=Gene)) + geom_line() + facet_wrap(~Chr,ncol=1)




gem <- cbind(c("chr10","chr10","chr2","chr2","chr2","chr2","chr3","chr3"),c(69905300,71098195,15883478,18112473,16745172,18767451,108565815,111263826),rep(1,8),rep('..',8),rep('Gem',8))
gem <- data.table(gem)
colnames(gem)<-c("Chr","Pos","HomozygousState","Quality","ROH")
gem$Chr <- factor(gem$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
gem$Pos <- as.numeric(gem$Pos)



# plotting with the new blockFinderAppend function
CCGO59bf <- blockFinderAppend(CCGO59,"59")
CCGO60bf <- blockFinderAppend(CCGO60,"60")
CCGO61bf <- blockFinderAppend(CCGO61,"61")
CCGO62bf <- blockFinderAppend(CCGO62,"62")



ggplot(data=CCGO59_f, aes(x=Pos,y=HomozygousState)) + geom_segment(data=subset(CCGO59_f,Class=="Block"),aes(x=Pos,xend=End,y=1,yend=1),colour="Red",size=5)  + geom_line() + facet_wrap(~Chr,ncol=1)

test <- rbind(CCGO59bf,CCGO60bf,CCGO61bf,CCGO62bf)
ggplot(data=test, aes(x=Pos,y=HomozygousState,colour=Subject)) + geom_segment(data=subset(test,Class=="Block"&Subject==59),aes(x=Pos,xend=End,y=1,yend=1),colour="Red",size=5) +geom_segment(data=subset(test,Class=="Block"&Subject==62),aes(x=Pos,xend=End,y=0,yend=0),colour="Blue",size=5)  + geom_line() + facet_wrap(~Chr,ncol=1) + geom_point()
ggplot(data=test, aes(x=Pos,y=HomozygousState,colour=Subject)) + 
  geom_segment(data=subset(test,Class=="Block"&Subject==59),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Red",size=5) +
  geom_segment(data=subset(test,Class=="Block"&Subject==60),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Green",size=5) +
  geom_segment(data=subset(test,Class=="Block"&Subject==61),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Orange",size=5) +
  geom_segment(data=subset(test,Class=="Block"&Subject==62),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Blue",size=5) + 
  facet_wrap(~Chr,ncol=2)