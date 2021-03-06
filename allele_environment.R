#!/usr/bin/env Rscript

# Script to take in chr-split VCF files (tested only on ExAC 0.3.1)
# and spit out the ExAC vcf with alleles, AFs, and distances for
# each variant (within a 500bp window)

# First keep only PASS variants (gzcat for mac, zcat for linux, PC? uh, I don't know)
# gzcat ~/Desktop/exac/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz | grep 'PASS\|^#' > ~/Desktop/exac/ExAC_nonTCGA.r0.3.1.sites.vep.PASSonly.vcf
# splits ExAC vcf by chr
# gzcat ExAC_nonTCGA.r0.3.1.sites.vep.PASSonly.vcf.gz | awk '!/^#/{print>>$1;close($1)}'

library(data.table)
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=1) {
  stop("Only one argument can be supplied (input file).\n", call.=FALSE)
}
# read in vcf which has already been split by chr and is bereft of the header
chr_vcf <- fread(args[1])
chr_vcf<-data.frame(chr_vcf)
colnames(chr_vcf) <- c("Chr","Position","ID","Ref","Alt","Qual","Filter","Info")
chr_vcf <- chr_vcf[!duplicated(chr_vcf),]
# pull AF into a new column
chr_vcf$AF <- apply(chr_vcf,1,function(x) {
  AF_pos <- grep('AF=',strsplit(x[8],';')[[1]])
  strsplit(x[8],";")[[1]][AF_pos]
})
# remove AF=
chr_vcf$AF <- sapply(chr_vcf$AF, function(x) substr(x,4,4000))
chr_vcf <- chr_vcf[,c(1,2,3,4,5,6,7,9)]
chr_vcf <- data.table(chr_vcf)  # may speed things up a bit, being a data.table
positions <- chr_vcf[,Position]
# chr_vcf[Position<(239313+250) & Position>(239313-250),]

# Two series of loops, the first finding all indices (genomic positions)
# that are within 250bp upstream of the variant. The second goes 250bp
# the other way. Very fast (~25s on chr3 on my i7 SSD 2013 15'' MacBook Pro)
#
# pre-initialize vector to save a bit of time
indices_up <- vector(mode='integer',length=length(positions))
system.time(
for (i in 1:length(positions)) {
  counter <- 1
  while ((i+counter)<length(positions) & 
         abs(positions[i]-positions[i+counter])<250) {
    counter=counter+1
  }
  indices_up[i] <- counter - 1
})
print("Finished with up")
indices_down <- vector(mode='integer',length=length(positions))
system.time(
for (i in length(positions):1) {
  counter <- 0
  while (
    abs(positions[i]-positions[i-counter])<250 ) {
    counter=counter+1
    if (i-counter<1) break #kills while loop when start position is going 0 or neg
  }
  indices_down[i] <- counter - 1
}) 
print("Finished with down")

# glue up and down together
indices <- cbind(indices_up,indices_down)

# now the slow part (~1hr on chr3)
# Using the indices above, Get all of the alleles
# Allele Frequencies (AF) and distances from the 
# variant
Alt_Alleles<-vector(mode='integer',length=nrow(chr_vcf))
Alt_Alleles<-list(Alt_Alleles)
AFs<-vector(mode='integer',length=nrow(chr_vcf))
AFs<-list(AFs)
Distances<-vector(mode='integer',length=nrow(chr_vcf))
Distances<-list(Distances)
system.time(
for (i in 1:nrow(chr_vcf)){
    start = i-indices[i,2]
    stop = i+indices[i,1]
    Alt_Alleles[i] <- list(unlist(strsplit(chr_vcf[start:stop,Alt],",")))
    AFs[i] <- list(as.numeric(unlist(strsplit(chr_vcf[start:stop,AF],","))))
    Distances[i] <- list(chr_vcf[start:stop,Position] - chr_vcf[i,Position])
} )

# first save to R data type
temp <- cbind(Alt_Alleles,AFs,Distances)
chr_vcf_R <- cbind(chr_vcf,temp)
output_name <- paste("VariantEnvironment_",args[1],sep='')
output_name_R <- paste(output_name,".Rdata",sep='')
save(chr_vcf_R,file=output_name_R)

# now collapse lists to save to plaintext file. DONT BOTHER. HAVE TO RECOMBINE ANY WAY. EASIER WITH RDATA
#Alt_Alleles<-sapply(Alt_Alleles,function(x) paste(x,collapse=','))
#AFs <- sapply(AFs,function(x) paste(x,collapse=',')) 
#Distances<-sapply(Distances,function(x) paste(x,collapse=','))

#temp <- cbind(Alt_Alleles,AFs,Distances)
#chr_vcf_P <- cbind(chr_vcf,temp)
#output_name_P <- paste(output_name,".txt",sep='')
#write.table(chr_vcf_P,file=output_name_P,quote=FALSE,sep='\t',row.names=FALSE)
