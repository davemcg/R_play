library(data.table)
# gzcat ExAC_nonTCGA.r0.3.1.sites.vep.PASSonly.vcf.gz | awk '!/^#/{print>>$1;close($1)}'
# splits ExAC vcf by chr
chr3 <- fread('~/Desktop/exac/3')
chr3<-data.frame(chr3)
colnames(chr3) <- c("Chr","Position","ID","Ref","Alt","Qual","Filter","Info")
chr3$AF <- apply(chr3,1,function(x) strsplit(x[8],";")[[1]][12])
chr3 <- chr3[,c(1,2,3,4,5,6,7,9)]
