# grab output files from allele_environment.R and merge into on file

#library(data.table)
#library(dplyr)

setwd('/Volumes/mcgaugheyd/exac/') # mount cyclops home dir
#files to slurp in
input_files <- list.files(pattern='*.Rdata')

load(input_files[length(input_files)])
temp <- head(chr_vcf_R,1)
for (i in input_files){
  load(i)
  temp<- rbind(temp,chr_vcf_R)
}

exac0.3.1_variantEnvironment <- temp[-1,]

# Quick and dirty grab the first AF (which is the most common) for each variant line
exac0.3.1_variantEnvironment$AFn <- sapply(exac0.3.1_variantEnvironment$AF,function(x) strsplit(x,',')[[1]][1])
exac0.3.1_variantEnvironment$AFn <- as.numeric(exac0.3.1_variantEnvironment$AFn)

library(ggplot2)
library(data.table)
library(dplyr)

#> summary(exac0.3.1_variantEnvironment$AFn)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.0000094 0.0000094 0.0000188 0.0045970 0.0000394 1.0000000 
### Odd, shouldn't have any AF above 1

# > table(exac0.3.1_variantEnvironment$AFn>0.5)
# FALSE    TRUE 
# 7630556   24623 
### 24,623 out of 7.6 million. Positions where the Reference should be updated?

# I'll just take the rest for now

exac <- dplyr::filter(exac0.3.1_variantEnvironment,AFn < 0.5)

exac <- exac %>% mutate(quantile = ntile(AFn, 10)) # get quantiles for AF
# now need to expand the distances



exac$Num_Surrounding <- sapply(exac$Distances,function(x) length(x))
quantiles <- unlist(apply(exac,1,function(x) rep(x['quantile'],x['Num_Surrounding'])))
distances <- unlist(exac$Distances)

dq <- data.table(cbind(distances,quantiles))
dqS <- dq %>% sample_n(5000000)
ggplot(data=dqS,aes(x=distances,colour=factor(quantiles))) + geom_density()
ggplot(data=dqS,aes(x=distances,colour=factor(quantiles))) + geom_freqpoly(binwidth=1) + scale_colour_brewer(palette="PRGn") + theme_dark() + scale_y_continuous(limits=c(0,30000))
