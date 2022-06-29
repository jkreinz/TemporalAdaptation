#!/usr/bin/env Rscript
itt = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(data.table)

cg_012<-fread(paste("../cg_ordered_",itt,"_matches.012",sep=""),header=F,sep=" ",strip.white = TRUE)
cg_012$V1 <- gsub("\\[b37\\].*","",cg_012$V1)
cg_names<-fread("cg.ind")
names(cg_012)<-c("snp",cg_names$V1)

herb_012<-fread(paste("../herb_ordered_",itt,"_matches.012",sep=""),header=F,sep=" ",strip.white = TRUE)
herb_012$V1<- gsub("\\[b37\\].*","",herb_012$V1)
herb_names<-fread("herb.ind",header=F)
names(herb_012)<-c("snp",herb_names$V1)

herb_cg_matched<- inner_join(herb_012,cg_012,by="snp") %>% select(-snp)
herb_cg_matched<- as.data.frame(t(herb_cg_matched))
herb_cg_matched$sample<-row.names(herb_cg_matched)

info<-read.table("3waymerged_sampleinfo.txt",sep="\t",header = T)

herb_metadata_012<-inner_join(herb_cg_matched,info,by="sample")
long_012<-melt(herb_metadata_012,id.vars=c("sample","env","sex","lat","long","year","state"))

pval<-list()
zval<-list()
slope<-list()
afmin<-list()
afmax<-list()
slope_error<-list()

ag<-herb_metadata_012 %>% filter(env!="Nat") %>% filter(year > 1869)
newdata <- data.frame(year = c(1870, 2018))
year=as.numeric(ag$year)
for (i in 1:154) {
  geno<-as.numeric(ag[[i]])/2
  lm1<-glm(geno ~ year, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)

  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  slope_error[[i]]<-test$`Std. Error`[2]

  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1]
  afmax[[i]]<-probabilities[2]

}

results_ag<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error),afmin=unlist(afmin),afmax=unlist(afmax))
mean(results_ag$afmax)
results_ag$variable<-paste("snp",seq(1:154),sep = "_")
long_012<-long_012[long_012$year > 1869,]
long_012$value<-long_012$value/2

fulltime_ag<-glm(data=long_012, value ~ year + variable, family="binomial")
#summary(fulltime_ag)

pre_ag<-glm(data=long_012[long_012$year < 1960,], value ~ year + variable, family="binomial")
post_ag<-glm(data=long_012[long_012$year > 1959,], value ~ year + variable, family="binomial")

percentin<-sum(results_ag$slope>0)/154
meanchange<-mean((results_ag$afmax-results_ag$afmin))
joint_s<-fulltime_ag$coefficients[2]
joint_s_pre<-pre_ag$coefficients[2]
joint_s_post<-post_ag$coefficients[2]
medianchange<-median((results_ag$afmax-results_ag$afmin))
meanslope<-mean(results_ag$slope)
medianslope<-median(results_ag$slope)
meansloperr<-mean(results_ag$slope_error)
mediansloperr<-median(results_ag$slope_error)
meanpval<-mean(results_ag$pval)
medianpval<-median(results_ag$pval)

write.table(t(c(percentin, meanchange, joint_s, joint_s_pre, joint_s_post, medianchange, meanslope, medianslope, meansloperr, mediansloperr, meanpval, medianpval)), paste("summarystat_",itt,".txt",sep=""),quote=F,col.names=F, row.names=F)
