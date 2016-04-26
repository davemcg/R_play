totalAZ <- function(the_df,num_of_var,genomicDelta) {
  colnames(the_df) <- c("Chr","Pos","HomozygousState","Quality")
  # select chroms to loop through
  chroms <- unique(the_df$Chr)
  # setup block df for ROH blocks with dummy values to be removed later
  runs <- data.table(rbind(c(0,0,0,0,0,0)))
  colnames(runs) <- c('values','starts','stops','lengths','GenomicDistance','Chr')
  for (i in chroms){
    chr_df <- subset(the_df,Chr=="chr1")
    # rle2 will ID consecutives runs and give the indices for the runs
    chr_runs <- data.table(rle2(chr_df$HomozygousState,indices=TRUE))
    chr_runs$GenomicDistance<-apply(chr_runs,1,function(x) abs(chr_df[x[2]]$Pos - chr_df[x[3]]$Pos))
    chr_runs$Chr <- i
    runs <- rbind(runs,chr_runs)
  }
  runs <- runs[-1,]
  total_autozygous_length <- data.frame(runs %>% group_by(values) %>% summarise(sum=sum(GenomicDistance)))[2,2]
  total_hw_length <- data.frame(runs %>% group_by(values) %>% summarise(sum=sum(GenomicDistance)))[1,2]
  output <- c(total_autozygous_length,total_hw_length)
  return(output)
}



CCGO59 <- fread('~/Dropbox/CCGO_800059.bcftools.rohViterbit',skip=3)
CCGO60 <- fread('~/Dropbox/CCGO_800060.bcftools.rohViterbi',skip=3)
CCGO61 <- fread('~/Dropbox/CCGO_800061.bcftools.rohViterbi',skip=3)
CCGO62 <- fread('~/Dropbox/CCGO_800062.bcftools.rohViterbi',skip=3)

CCGO159v <- fread('~/Dropbox/CCGO_800159.bcftools.rohViterbi',skip=3)
CCGO160v <- fread('~/Dropbox/CCGO_800160.bcftools.rohViterbi',skip=3)
CCGO161v <- fread('~/Dropbox/CCGO_800161.bcftools.rohViterbi',skip=3)