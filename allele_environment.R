library(data.table)
# gzcat ExAC_nonTCGA.r0.3.1.sites.vep.PASSonly.vcf.gz | awk '!/^#/{print>>$1;close($1)}'
# splits ExAC vcf by chr
chr3 <- fread('~/Desktop/exac/3')
chr3<-data.frame(chr3)
colnames(chr3) <- c("Chr","Position","ID","Ref","Alt","Qual","Filter","Info")
chr3$AF <- apply(chr3,1,function(x) strsplit(x[8],";")[[1]][12])
chr3 <- chr3[,c(1,2,3,4,5,6,7,9)]
chr3 <- chr3[!duplicated(chr3),]
chr3 <- data.table(chr3)
# chr3[Position<(239313+250) & Position>(239313-250),]

countAdjacentAlleles <- function(chr_data_table, window_size) {
  # too damn slow. About 8s to do the first 10,000 with 250bp window in chr3
  #chr_data_table <- data.table(chr_data_table)
  out <- sapply(chr_data_table[,Position],function(x) 
    nrow( 
      chr_data_table[Position < (x + window_size) & 
                       Position > (x - window_size),])
  )
  return(out)
}

countAdjacentAlleles2 <- function(chr_vcf, window_size) {
  # a bit faster
  chr_vcf <- data.frame(chr_vcf)
  out <- sapply(chr_vcf[,'Position'],function(x) 
    nrow( 
      subset(chr_vcf, Position < (x + window_size) & 
                       Position > (x - window_size))
      )
  )
  return(out)
}

countAdjacentAlleles3 <- function(chr_vcf, window_size) {
  # a bit faster
  chr_vcf <- data.table(chr_vcf)
  positions <- as.integer(chr_vcf$Position)
  out <- vector(mode="integer",length=nrow(chr_vcf))
  out <- sapply(positions,function(x) 
    nrow( 
      chr_vcf %>% filter(Position < (x + window_size) & 
                           Position > (x - window_size))
    )
  )
  return(out)
}

total <- vector(mode='integer',length=length(positions))
system.time(
for (i in 1:length(positions)) {
  counter <- 1
  while ((i+counter)<length(positions) & 
         abs(positions[i]-positions[i+counter])<250) {
    counter=counter+1
  }
  total[i] <- counter - 1
} )

total_down <- vector(mode='integer',length=length(positions))
system.time(
for (i in length(positions):1) {
  counter <- 0
  while (
    abs(positions[i]-positions[i-counter])<250 ) {
    counter=counter+1
    if (i-counter<1) break #kills while loop when hitting the end
  }
  total_down[i] <- counter - 1
} )

alleles<-vector(mode='integer',length=nrow(chr3))
alleles<-list(alleles)
AFs<-vector(mode='integer',length=nrow(chr3))
AFs<-list(AFs)
distances<-vector(mode='integer',length=nrow(chr3))
distances<-list(distances)
system.time(
for (i in 1:nrow(chr3)){
    start = i-both[i,2]
    stop = i+both[i,1]
    alleles[i] <- list(unlist(strsplit(chr3[start:stop,Alt],",")))
    AFs[i] <- list(as.numeric(unlist(strsplit(chr3[start:stop,AF],","))))
    distances[i] <- list(chr3[start:stop,Position] - chr3[i,Position])
} )
  
  