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
pdf('CCGO_05-16-16_relatedness.pdf',width=8,height=6)
ggdendrogram(hc, size=12)
dev.off()


# color labelling fams
# http://stackoverflow.com/questions/8045538/labelling-ggdendro-leaves-in-multiple-colors
x <- as.matrix(rel)
dd.row <- as.dendrogram(hclust(dist(t(x))))
ddata_x <- dendro_data(dd.row)
p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
labs <- label(ddata_x)
labs <- labs %>% arrange(as.character(label))
labs$Family <- c(rep("CCGO_FAM_800016and18",5),rep("CCGO_FAM_800044",3),rep("CCGO_FAM_800062",4),rep("CCGO_FAM_800067and71",5),rep("CCGO_FAM_800118",3),rep("CCGO_FAM_800157",4),rep("CCGO_FAM_800160",3),rep("CCGO_FAM_800298",4),"CCGO_FAM_800308","CCGO_FAM_800347")
labs <- labs %>% arrange(x)
p2 + geom_text(data=label(ddata_x),
               aes(label=label, x=x, y=-0.1, colour=labs$Family, angle=90)) + scale_y_continuous(expand=c(0.05,0.05)) +
               theme(axis.ticks.x = element_blank(),axis.text.x=element_blank()) + xlab("") + ylab("Relatedness (0 is identical)")
