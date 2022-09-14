#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
i=args[1]

library(data.table)

genos<-read.table("m2_allsamples_nomulti.012")
inds<-read.table("m2_allsamples_nomulti.012.indv")
snps<-read.table("m2_allsamples_nomulti.012.pos")

year<-gsub("i","",gsub("t","",(gsub("_sample_i0","",inds$V1))))

genos_long<-melt(genos[,-1])

year2<-as.numeric(unlist(lapply(strsplit(year,"_"), `[[`, 1)))
genos_long$year<-rep(year2,ncol(genos)-1)

library(dplyr)
genos_long$value<-genos_long$value/2

#lm1<-genos_long %>% glm(formula=value~year+variable, family="binomial")
#summary(lm1)

newdata <- data.frame(year = c(0, 141))
pval<-list()
zval<-list()
slope<-list()
afmin<-list()
afmax<-list()
slope_error<-list()

for (i in 2:ncol(genos)) {
  geno<-genos[[i]]/2
  lm1<-glm(geno ~ year2, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)

  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  slope_error[[i]]<-test$`Std. Error`[2]

  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1]
  afmax[[i]]<-probabilities[2]

}

results_sims<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error),afmin=unlist(afmin),afmax=unlist(afmax))
selcoef<-read.table("selection_coefficients.txt")
results_sims$real_s<-selcoef$V1
results_sims$slope_2<-(results_sims$slope*2)
results_sims$slope_error_2<-(results_sims$slope_error*2)
results_sims_dropouts<-results_sims %>% filter(slope<5 & slope >-1) #removing true outliers
cor(results_sims_dropouts$slope_2,results_sims_dropouts$real_s)
