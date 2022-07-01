#analyses of finescale ancestry across the genome, from LAMP

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

anc<-fread("ancwindows_byscaf_bysample_20snpthresh.txt") #contemp windowed filt

#metadata
sampinfo<-fread("contemp_wpair.txt")
sampinfo<-sampinfo[,c(1,2,4)]
names(sampinfo)<-c("variable","env","pair")
sampinfo$variable<-as.character(as.integer(sampinfo$variable))

coord<-fread("contemp_coords.txt")
names(coord)<-c("variable","long","lat")
coord$variable<-as.character(as.integer(coord$variable))
moreinfo<-inner_join(coord,sampinfo,by="variable")

#merge
anc<-inner_join(anc,moreinfo,by="variable")

############
#Figure S10
############
#diff in ancestry by env across the genome

anc$scaf<-gsub("scaf","",anc$scaf)
anc$scaf<-as.numeric(anc$scaf)
anc %>% group_by(win, env, scaf) %>% 
  dplyr::summarise(mean_anc=mean(value), anc_std=stde(value)) %>% 
  ggplot(aes(win/10, mean_anc, color=env)) +
  geom_point(size=.6) +
  geom_path() +
  ylim(0.5,1) +
  ylab("var. rudis ancestry proportion") +
  theme_bw() +
  #geom_vline(xintercept = 22.5, lty="dashed") +
  #geom_vline(xintercept = 30, lty="dashed") +
  xlab("Genomic Position (Mb)") +
  geom_errorbar(aes(ymin=mean_anc-anc_std, ymax=mean_anc+anc_std),alpha=.5) +
  facet_wrap(~scaf, scales="free_x") 


############
#Figure 4D
############

#merge ancestry with XPEHH, gene density, and recombination rate info

byenv1 <- anc %>% group_by(win, scaf, pair, env) %>% dplyr::summarise(mean_anc=mean(value), long = first(long))
#get diff in ancestry by pair
ag <- byenv1[byenv1$env == "Ag",]
nat <- byenv1[byenv1$env == "Nat",]
ag_nat<-inner_join(ag, nat, by=c("win","scaf","pair"))
ag_nat$anc_diff<- ag_nat$mean_anc.x - ag_nat$mean_anc.y

test<-melt(ag_nat[,c(1,2,3,4,10)], id.vars=c("win","scaf","pair","anc_diff"))


#recombination rate data
recomb_winds<-fread("100kb_windowed_recomb.txt")
names(recomb_winds)<-c("win","median_recomb","mean_recomb","n_snps","start","end","scaf")
recomb_winds$scaf<-gsub(pattern = "Scaffold_",replacement = "",recomb_winds$scaf)
recomb_winds$scaf<-as.numeric(recomb_winds$scaf)
env_diff_bypair$scaf<-as.numeric(gsub(pattern = "scaf",replacement = "",env_diff_bypair$scaf))
head(recomb_winds)
head(env_diff_bypair)

ag_nat<-inner_join(ag, nat, by=c("win","scaf","pair"))
#ag_nat$anc_diff<- ag_nat$mean_anc.x - ag_nat$mean_anc.y
test<-melt(ag_nat[,c(1,2,3,5,8)], id.vars=c("win","scaf","pair"))
head(test)

#get the mean pairwise diff in ancestry across envs
anc_byenv<-test %>% group_by(win, scaf,variable) %>% dplyr:::summarise(anc=mean(value,na.rm=T))
head(anc_byenv)
anc_byenv$scaf<-as.numeric(gsub(pattern = "scaf",replacement = "",anc_byenv$scaf))

#merge mean diff pairwise ancestry and recombination rate in windows
anc_recomb<-inner_join(anc_byenv,recomb_winds, by=c("scaf","win"))
head(anc_recomb)

#import gene density
gene_dens<-fread("genomewide_genedensities.txt")
names(gene_dens)<-c("scaf","win","end","genes_n","gene_prop")  
gene_dens$scaf<-as.numeric(gsub("Scaffold_","",gene_dens$scaf))
gene_dens_bigger<- gene_dens %>% mutate(win=floor((win)/100000)) %>%
  group_by(win, scaf) %>% 
  dplyr::summarise(mean_dens=mean(gene_prop), median_dens=median(gene_prop),gene_n=sum(genes_n))

#merge gene density with anc and recomb
anc_recomb$scaf<-as.numeric(as.character(anc_recomb$scaf))
anc_recomb_dens<-inner_join(anc_recomb,gene_dens_bigger, by=c("scaf","win"))
anc_recomb_dens <- anc_recomb_dens %>% filter(mean_recomb < 50) #remove outliers

#import CMH
cmh<-fread("~/ag_v_nat_withinpairs_maf01.cmh")
cmh_wind<- cmh %>% mutate(win=floor((BP)/100000)) %>%
  group_by(win, CHR) %>% 
  dplyr::summarise(mean_cmh=mean(P,na.rm = T), median_cmh=median(P,na.rm = T),cmh_n=n())

cmh_wind$scaf<-as.numeric(gsub("Scaffold_","",cmh_wind$CHR))
head(cmh_wind)
head(anc_recomb_dens)
#merge CMH with anc, recomb, and dens
cmh_anc_recomb_dens<-inner_join(cmh_wind, anc_recomb_dens, by=c("scaf","win") )

#import XPEHH
xpehh2<-fread("XPEHH_nat_ag_allchroms.txt")
head(xpehh2)
names(xpehh2)[c(1,2,10)]<-c("scaf","pos","xpehh")
xpehh2$scaf<-as.numeric(xpehh2$scaf)
xpehh_wind<- xpehh2 %>% mutate(win=floor((as.numeric(pos))/100000)) %>%
  group_by(win, scaf) %>% 
  dplyr::summarise(mean_xpehh=mean(xpehh,na.rm = T), median_xpehh=median(xpehh,na.rm = T),xpehh_n=n())

#merge XPEHH with gene density, recomb, and anc
anc_fst_dens_recomb_xpehh2<-inner_join(xpehh_wind, cmh_anc_recomb_dens, by=c("scaf","win") )
head(anc_fst_dens_recomb_xpehh2)

#run multiple regression
lm3<-lm(data=anc_fst_dens_recomb_xpehh2, anc ~  scaf + mean_dens + mean_xpehh * variable + mean_cmh * variable + mean_recomb + variable ) #circular to compare envfst to diff in env anc
summary(lm3)
Anova(lm3,type=3)

#check least squares mean difference in ancestry across environments
library(lsmeans)
lsmeans(lm3, "variable") #

#finally, Figure 4D

xpehh_plot<-ggplot(anc_fst_dens_recomb_xpehh2, aes(mean_xpehh,anc, color=variable)) +
  geom_point(alpha=.5) +
  geom_smooth(method="lm") +
  theme_bw() +
  scale_color_manual(labels=c("Ag","Nat"),values=c("grey","grey40")) +
  labs(y="Proportion of var. rudis ancestry", x="Ag-Nat XPEHH, 100kb average",
       color="Environment") 

cmh_plot<-ggplot(anc_fst_dens_recomb_xpehh2, aes(-log10(mean_cmh),anc, color=variable)) +
  geom_point(alpha=.5) +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(y="Proportion of var. rudis ancestry", x="-log10(CMH p value), 100 kb average",
       color="Environment") +
  scale_color_manual(labels=c("Ag","Nat"),values=c("grey","grey40"))

library(cowplot)
plot_grid(cmh_plot, xpehh_plot, ncol=2,rel_widths = c(.5,.5))

