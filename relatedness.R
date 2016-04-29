#/bin/bash/Rscript

# plots dendrogram of relatedness between samples in a multi-sample VCF
# out.relatedness2 comes from vcftools --relatendess2
# usage: vcftools --gzvcf CCGO.b37.bwa-mem.hardFilterSNP-INDEL.vcf.gz.VEP.GRCh37.vcf.gz --relatedness2

library(data.table)
library(reshape2)
library(ggdendro)

relations <- data.frame(fread('~/Desktop/out.relatedness2'))
relate <- relations[,c("INDV1","INDV2","RELATEDNESS_PHI")]
rel <- dcast(relate, INDV1~INDV2)
row.names(rel)<-rel$INDV1
rel <- rel[,-1]

d <- dist(as.matrix(rel))
hc <- hclust(d)
pdf('test.dpf',width=8,height=6)
ggdendrogram(hc, size=12)
dev.off()