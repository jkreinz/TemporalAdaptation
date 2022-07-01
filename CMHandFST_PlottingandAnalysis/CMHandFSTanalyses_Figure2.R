cmh<-fread("~/ag_v_nat_withinpairs_maf01.cmh")
cmh$fdr <- p.adjust(cmh$P,method="fdr")
cmh$bon <- p.adjust(cmh$P,method="bonferroni")
fdr_sig <- cmh[cmh$fdr < 0.1,]
bon_sig <- cmh[cmh$bon < 0.1,]

#clumped data with focal SNP indicated
clump<-fread("ag_v_nat_withinpairs_maf01_1mb_r25.clumped")
cmh2<-cmh[,c(1,3,8)]
clump2<-clump[,c(1,4,5)]
cmh2$type<-"inlier"
clump2$type<-"indexSNP"
clump2[clump2$P < 5.55e-06,]$type<-"fdrsig_indexSNP"

#merge genome wide CMH with clumped CMH
cmh_clump<-rbind(cmh2,clump2)
cmh_clump$type<-as.factor(cmh_clump$type)
levels(cmh_clump$type)
cmh_clump$CHR<-gsub("Scaffold_", "", cmh_clump$CHR)
cmh_clump$CHR<- as.numeric(cmh_clump$CHR)

table(cmh_clump$type)
#1) number of LD thinned, FDR sig snps, 
#2) number of index (clumped, LD thinned peaks)
#3) all other SNPs across the genome

###########
#figure 2B
###########

cmh_clump %>% filter(CHR=="11") %>% 
  ggplot(aes(BP/1000000,-log10(P), fill=type, color=type, alpha=type)) +
  # facet_grid(. ~ CHR, scales = "free_x", space = "free_x") +
  geom_vline(xintercept = 23.999751,lwd=.75) +
  geom_vline(xintercept = 24.250051,col="grey45",lwd=.5) +
  geom_point() +
  theme_bw() +
  xlab("Genomic Position") +
  geom_hline(yintercept = -log10(max(fdr_sig$P)),lty="dashed") + 
  geom_hline(yintercept = -log10(max(bon_sig$P)),lty="dashed",color="lightgrey") + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("CMH p-value") +
  scale_color_manual(values=c("firebrick","#421513","black")) +
  scale_fill_manual(values=c("firebrick","#421513","black")) +
  scale_alpha_manual(values=c(1,.2,.2)) +  
  theme(legend.position = "none") +
  xlim(20,30)

##############
#figure 2A
##############

cmh_plot<-cmh_clump %>% #filter(CHR=="Scaffold_11") %>% 
  ggplot(aes(BP/1000000,-log10(P), fill=type, color=type, alpha=type)) +
  facet_grid(. ~ CHR, scales = "free_x", space = "free_x") +
  geom_point() +
  theme_bw() +
  xlab("Genomic Position") +
  geom_hline(yintercept = -log10(max(fdr_sig$P)),lty="dashed") + 
  geom_hline(yintercept = -log10(max(bon_sig$P)),lty="dashed",color="lightgrey") + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("CMH p-value") +
  scale_color_manual(values=c("firebrick","#421513","black")) +
  scale_fill_manual(values=c("firebrick","#421513","black")) +
  scale_alpha_manual(values=c(1,.2,.2)) +
  theme(legend.position = "none") 


#############
#figure 2C
#############

library(data.table)
fst<-fread("ag_v_nat.weir.fst",header = T,na.strings = c("NaN","-Inf"))
fst[fst$WEIR_AND_COCKERHAM_FST < 0]<-0

cands<-read.table("Desktop/herb AF plots/canadidate_fst.txt")

ggplot() +
  geom_histogram(data=fst,aes(WEIR_AND_COCKERHAM_FST), fill="lightgrey",color="black") +
  theme_classic() +
  labs(x=expression(paste(F[ST])),y="Count") +
  geom_vline(xintercept = cands$V3,lty="dashed")



##############
#exploring how CMH correlates with FST
##############

#match fst (no filt) and CMH (maf > 0.01)
fst<-fread("ag_v_nat.weir.fst",header = T,na.strings = c("NaN","-Inf"))
names(fst)<-c("CHR","BP","P")
fst$CHR<-as.numeric(gsub("Scaffold_","",fst$CHR))
cmh$CHR<-as.numeric(gsub("Scaffold_","",cmh$CHR))
fst_less<-inner_join(fst, cmh, by=c("CHR","BP"))

#again, import clumped results to identify focal SNPs
clump<-fread("ag_v_nat_withinpairs_maf01_1mb_r25.clumped")
clump2<-clump[,c(1,4,5)]
fst_less$type<-"inlier"
clump2$type<-"indexSNP"
clump2[clump2$P < 5.55e-06,]$type<-"fdrsig_indexSNP" #fdr p-value thresh
fst_less<-fst_less[,c(1,2,3,4)]
names(fst_less)[3]<-"P"

clump2[clump2$type == "fdrsig_indexSNP",]
clump2$CHR<-gsub("Scaffold_","",clump2$CHR)
clump2$CHR<-as.numeric(clump2$CHR)
names(fst_less)[3]<-"fst"
fst_clump<-left_join(fst_less, clump2, by=c("CHR","BP"))
head(fst_clump)

fst_clump$type[is.na(fst_clump$P)] <- "inlier"
fst_clump$type<-as.factor(fst_clump$type)
levels(fst_clump$type)
table(fst_clump$type)

fst_plot<-ggplot(data=fst_clump, aes(-log10(P),Fst, color=type, fill=type, alpha=type)) +
  geom_point() +
  scale_color_manual(values=c("firebrick","black","black")) +
  scale_fill_manual(values=c("firebrick","black","black")) +
  scale_alpha_manual(values=c(1,.2,.2)) +
  geom_hline(yintercept = quantile(fst_clump$fst, probs=0.999),lty="dashed") + 
  geom_vline(xintercept = -log10(max(fdr_sig$P)),lty="dashed") +   
  theme_bw() +
  theme(legend.position = "none") +
  labs(x="-log10(CMH p-value)")

cowplot::plot_grid(cmh_plot, fst_plot, align = T, ncol = 1, nrow=2)


#covariance in focal SNPs?
pos1<-read.table("herb_cmh_FDR10p_clumped_matched_refalt_matched.012.pos")
pos2<-read.table("herb_cmh_FDR10p_clumped_matched_refalt_toflip.012.pos")
fullpos<-rbind(pos1,pos2)
fullpos$V1<-as.numeric(gsub("Scaffold_","",fullpos$V1))
names(fullpos)<-c("CHR","BP")
head(fullpos)

fst<-fread("ag_v_nat.weir.fst",header = T,na.strings = c("NaN","-Inf"))
mean_fst<-fst
head(mean_fst)
names(mean_fst)<-c("CHR","BP","fst")
hist(mean_fst$fst)
head(cmh)
library(dplyr)

cmh$fdr <- p.adjust(cmh$P,method="fdr")
fdr_sig <- cmh[cmh$fdr < 0.1,]
mean_fst$CHR<-as.numeric(gsub("Scaffold_","",mean_fst$CHR))
cmh_pairedfst<-inner_join(cmh,mean_fst,by=c("CHR","BP"))

cmh_pairedfst$fst_top<- cmh_pairedfst$fst > quantile(cmh_pairedfst$fst,probs = 0.999)
cmh_pairedfst$fdr <- FALSE 
cmh_pairedfst$fdr[cmh_pairedfst$P < 5.55e-06] <- TRUE
table(cmh_pairedfst$fst_top, cmh_pairedfst$fdr)
#using a 0.999 Fst threshold
#9/402, or 0.022 focal SNPs outliers for CMH but not Fst
#393/402, or 0.977 focal SNP outliers for both CMH and Fst

cmh_pairedfst$fst_toptop<- cmh_pairedfst$fst > quantile(cmh_pairedfst$fst,probs = 0.9999)
#using a 0.9999 Fst threshold
table(cmh_pairedfst$fst_toptop, cmh_pairedfst$fdr)


#what is the correlation of CMH and Fst for all sites across the genome?
cor(-log10(cmh_pairedfst$P), cmh_pairedfst$fst, method = "spearman", use = "complete.obs")
cor(cmh_pairedfst$CHISQ, cmh_pairedfst$fst, method = "pearson", use = "complete.obs")

