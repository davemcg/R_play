library(accelerometry)
library(ggplot2)
library(data.table)
library(dplyr)

# read in output from bcftool roh
######## USE VITERBI training!!!############
CCGO59 <- fread('~/Dropbox/CCGO_800059.bcftools.roh2',skip=3)
CCGO60 <- fread('~/Dropbox/CCGO_800060.bcftools.roh2',skip=3)
CCGO61 <- fread('~/Dropbox/CCGO_800061.bcftools.roh2',skip=3)
CCGO62 <- fread('~/Dropbox/CCGO_800062.bcftools.roh2',skip=3)

CCGO159 <- fread('~/Dropbox/CCGO_800159.bcftools.roh',skip=3)
CCGO160 <- fread('~/Dropbox/CCGO_800160.bcftools.roh',skip=3)
CCGO161 <- fread('~/Dropbox/CCGO_800161.bcftools.roh',skip=3)

# this function use rle2 to identify consecutive positions with homozygosity
# and appends the blocks (currently coded for >100 variants and >1mb)
# to the inputted df
blockFinderAppend <- function(the_df, Subject,num_of_var,genomicDelta) {
  colnames(the_df) <- c("Chr","Pos","HomozygousState","Quality")
  the_df$End <- the_df$Pos 
  the_df$Subject <- Subject
  the_df$Class <- "Variant"
  the_df <- the_df %>% dplyr::select(Chr, Pos, End, HomozygousState, Quality, Subject, Class)
  # select chroms to loop through
  chroms <- unique(the_df$Chr)
  # setup block df for ROH blocks with dummy values to be removed later
  blocks <- data.table(rbind(c(0,0,0,0,0,0,0)))
  colnames(blocks) <- c('Chr', 'Pos', 'End', 'HomozygousState', 'Quality', 'Subject', 'Class')
  for (i in chroms){
    chr_df <- subset(the_df,Chr==i)
    # rle2 will ID consecutives runs and give the indices for the runs
    chr_runs <- data.table(rle2(chr_df$HomozygousState,indices=TRUE))
    # filter to keep runs with > 100 positions
    over <- chr_runs %>% dplyr::filter(((stops-starts) > num_of_var) & values==1) %>% arrange(-lengths)
    # calculate genomic space in each ROH
    over$GenomicDistance<-apply(over,1,function(x) abs(chr_df[x[2]]$Pos - chr_df[x[3]]$Pos))
    over$Pos <- apply(over,1,function(x) chr_df[x[2]]$Pos-0)
    over$End <- apply(over,1,function(x) chr_df[x[3]]$Pos-0)
    over$Chr <- i
    # only keep roh with >1,000,000kb
    over <- over %>% filter(GenomicDistance > genomicDelta)
    over$HomozygousState <- "1"
    over$Quality <- "0"
    over$Subject <- Subject
    over$Class <- "Block"
    over <- over %>% select(Chr, Pos, End, HomozygousState, Quality, Subject, Class)
    blocks <- rbind(blocks,over)
  }
  blocks <- blocks[-1,]
  output <- rbind(the_df, blocks)
  output$Pos <- as.numeric(output$Pos)
  output$End <- as.numeric(output$End)
  output$HomozygousState <- as.numeric(output$HomozygousState)
  # order chromosomes
  output$Chr <- factor(output$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
  
  return(output)
}

# plotting with the new blockFinderAppend function
CCGO59bf <- blockFinderAppend(CCGO59,"59")
CCGO60bf <- blockFinderAppend(CCGO60,"60")
CCGO61bf <- blockFinderAppend(CCGO61,"61")
CCGO62bf <- blockFinderAppend(CCGO62,"62")

CCGO159bf <- blockFinderAppend(CCGO159,"159")
CCGO160bf <- blockFinderAppend(CCGO160,"160")
CCGO161bf <- blockFinderAppend(CCGO161,"161")

ggplot(data=CCGO59bf, aes(x=Pos,y=HomozygousState)) + geom_segment(data=subset(CCGO59bf,Class=="Block"),aes(x=Pos,xend=End,y=1,yend=1),colour="Red",size=5)  + geom_line() + facet_wrap(~Chr,ncol=1)

test <- rbind(CCGO59bf,CCGO60bf,CCGO61bf,CCGO62bf,CCGO159bf,CCGO160bf,CCGO161bf)
test$Chr <- factor(test$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))

ggplot(data=test, aes(x=Pos,y=HomozygousState,colour=Subject)) + geom_segment(data=subset(test,Class=="Block"&Subject==59),aes(x=Pos,xend=End,y=1,yend=1),colour="Red",size=5) +geom_segment(data=subset(test,Class=="Block"&Subject==62),aes(x=Pos,xend=End,y=0,yend=0),colour="Blue",size=5)  + geom_line() + facet_wrap(~Chr,ncol=1) + geom_point()
ggplot(data=test, aes(x=Pos,y='Coverage',colour=Subject)) + 
  geom_point(size=0.1) + 
  geom_segment(data=subset(test,Class=="Block"&Subject==59),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Red",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==60),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Green",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==61),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Orange",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==62),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Blue",size=2) + 
  facet_wrap(~Chr,ncol=2) + theme_bw() + xlab('') + ylab('')

# hom rec variants for CCGO_800059-62
var <- cbind(c("chr12","chr7","chr1","chr15"),c(9465464,65548068,153314125,72190416),c(9465464,65548068,153314125,72190416),c(0,0,0,0),c('..','..','..','..'),c('62.var','62.var','62.var','62.var'),c('AR_variant','AR_variant','AR_variant','AR_variant'),c("RP11-22B23.1","ASL","PGLYRP4","MYO9A"))
colnames(var) <- c("Chr","Pos","End","HomozygousState","Quality","Subject","Class","Gene")
var <- data.table(var)
var$Chr <- factor(var$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
var$Pos <- as.numeric(var$Pos)
test2<-test
test2$Gene <- ''
test2 <- rbind(test2,var)
test2$End<-as.numeric(test2$End)
test2<-test2[!duplicated(test2),]
ggplot(data=test2, aes(x=Pos,y='Coverage',colour=Subject)) + 
  geom_point(size=0.1) + geom_point(data=subset(test2,Class=="AR_variant"),aes(x=Pos,y=Subject,shape=Gene)) +
  geom_segment(data=subset(test,Class=="Block"&Subject==59),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Red",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==60),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Green",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==61),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Orange",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==62),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Blue",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==159),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Yellow",size=2) + 
  geom_segment(data=subset(test,Class=="Block"&Subject==160),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Purple",size=2) +
  geom_segment(data=subset(test,Class=="Block"&Subject==161),aes(x=Pos,xend=End,y=Subject,yend=Subject),colour="Gray",size=2) +
  facet_wrap(~Chr,ncol=2) + theme_bw() + xlab('') + ylab('')




# pull in genotypes
gt <- fread("~/Dropbox/CCGO_800014-161.bwa-mem.hg19.GATK-3.4-46.raw.filterSNP-INDEL.copy.VEP.GRCh37.gt")

extractGT<- function(subjectVector){
  charGT <- sapply(subjectVector, function(x) strsplit(x,":")[[1]][1])
  # label unknown (contains '.') 
  charGT[grep('\\.',charGT)]<-NA
    #convert to numeric gt (0/0 = 0, 0/1 = 0.5, 1/1 = 1)
  charGT<-gsub('0\\/0',0,charGT)
  charGT<-gsub('0\\/1',0.5,charGT)
  charGT<-gsub('1\\/1',0,charGT)
  charGT<-as.numeric(charGT)
  return(charGT)
}

extractColN<- function(subjectVector,N){
  charV <- sapply(subjectVector, function(x) strsplit(x,":")[[1]][N])
  # label unknown (contains '.') 
  return(as.numeric(charV))
}

gt$CCGO59heterozygosity<-extractGT(gt$CCGO_800059)
colnames(heterozygosity)[1:2] <- c("Chr","Pos")
test59<-merge(CCGO59v_bf,gt,by=c("Chr","Pos"))
test59$Chr <- factor(test59$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))



CCGO59v <- fread('~/Dropbox/CCGO_800059.bcftools.rohViterbit',skip=3)
CCGO60v <- fread('~/Dropbox/CCGO_800060.bcftools.rohViterbi',skip=3)
CCGO61v <- fread('~/Dropbox/CCGO_800061.bcftools.rohViterbi',skip=3)
CCGO62v <- fread('~/Dropbox/CCGO_800062.bcftools.rohViterbi',skip=3)

CCGO59v_bf <- blockFinderAppend(CCGO59v,"59")
CCGO60v_bf <- blockFinderAppend(CCGO60v,"60")
CCGO61v_bf <- blockFinderAppend(CCGO61v,"61")
CCGO62v_bf <- blockFinderAppend(CCGO62v,"62")

test <- rbind(CCGO59v_bf,CCGO60v_bf,CCGO61v_bf,CCGO62v_bf)
test$Chr <- factor(test$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))










test59<-(blockFinderAppend(CCGO59v,'59',1,0))
test60<-(blockFinderAppend(CCGO60v,'60',1,0))
test61<-(blockFinderAppend(CCGO61v,'61',1,0))
test62<-(blockFinderAppend(CCGO62v,'62',1,0))

test <- rbind(test59,test60,test61,test62)
test$Chr <- factor(test$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))

var <- cbind(c("chr12","chr7","chr1","chr15"),c(9465464,65548068,153314125,72190416),c(9465464,65548068,153314125,72190416),c(0,0,0,0),c('..','..','..','..'),c('62.var','62.var','62.var','62.var'),c('AR_variant','AR_variant','AR_variant','AR_variant'),c("RP11-22B23.1","ASL","PGLYRP4","MYO9A"))
colnames(var) <- c("Chr","Pos","End","HomozygousState","Quality","Subject","Class","Gene")
var <- data.table(var)
var$Chr <- factor(var$Chr,levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"))
var$Pos <- as.numeric(var$Pos)
test2<-test
test2$Gene <- ''
test2 <- rbind(test2,var)
test2$End<-as.numeric(test2$End)
test2<-test2[!duplicated(test2),]




ggplot(data=test, aes(x=Pos,y=HomozygousState,colour=Subject)) + geom_line() + facet_grid(~Chr+Subject) + geom_point()

######################
# Viterbi data is WAY cleaner
#######################

ggplot(data=test2, aes(x=Pos,y='Coverage',colour=Subject)) + 
       geom_point(size=0.1) + geom_point(data=subset(test2,Class=="AR_variant"),aes(x=Pos,y=Subject,shape=Gene)) +
       geom_segment(data=subset(test,Class=="Block"&Subject=='59'),aes(x=Pos-500000,xend=End+500000,y=Subject,yend=Subject),colour="Red",size=2) +
       geom_segment(data=subset(test,Class=="Block"&Subject=='60'),aes(x=Pos-500000,xend=End+500000,y=Subject,yend=Subject),colour="Green",size=2) +
       geom_segment(data=subset(test,Class=="Block"&Subject=='61'),aes(x=Pos-500000,xend=End+500000,y=Subject,yend=Subject),colour="Orange",size=2) +
       geom_segment(data=subset(test,Class=="Block"&Subject=='62'),aes(x=Pos-500000,xend=End+500000,y=Subject,yend=Subject),colour="Blue",size=2) +
       facet_wrap(~Chr,ncol=2) + theme_bw() + xlab('') + ylab('')