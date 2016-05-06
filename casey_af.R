#
library(data.table)
vcf <- fread('~/Dropbox/second2.txt',skip=37)
# only keep info column with caseyTOPGG, caseySPPG, 1000G, ESP, ExAC AF
library(dplyr)
info <- vcf %>% select(INFO)

caseyTOPPG<-as.integer(sapply(info$INFO,function(x) (strsplit(strsplit(x,'\\|')[[1]][1],'=')[[1]][2])))
caseySPPG<-as.integer(sapply(info$INFO,function(x) strsplit(strsplit(strsplit(x,'\\|')[[1]][2],";")[[1]][1],'=')[[1]][2]))
gmaf <- as.numeric(sapply(info$INFO,function(x) (strsplit(strsplit(x,'\\|')[[1]][29],':')[[1]][2])))
exac <- as.numeric(sapply(info$INFO,function(x) (strsplit(strsplit(x,'\\|')[[1]][37],':')[[1]][2])))


gmaf[is.na(gmaf)] <- 1/1000000
exac[is.na(exac)] <- 1/1000000

caseyMAF <- caseyTOPPG / caseySPPG

AFs<-data.table(cbind(caseyTOPPG,caseySPPG,caseyMAF,gmaf,exac))

AFs2<-AFs[!is.na(caseyMAF)]
