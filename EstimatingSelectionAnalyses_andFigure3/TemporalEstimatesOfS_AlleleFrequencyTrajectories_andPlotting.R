#Estimate trajectory of allele frequencies from genotype data and infer selection, via logistic regression approach

#import genotypes for focal alleles in herbarium samples
library(data.table)
herb_012<-fread("~/herb_cmh_FDR10p_clumped_matched_refalt_flipped_final.012",na.strings = "-1") #output with plink
names(herb_012)<-paste("snp",seq(1:154),sep = "_")
herb_012<-herb_012/2 #rescale genotypes to allele frequencies i.e. 0,1,2, to 0,0.5,1
inds<-fread("~/herb_cmh_FDR10p_clumped_matched_refalt_matched.012.indv",header = F)
herb_012$sample<-gsub("HBO","HB0",inds$V1)

#import genotypes for focal alleles in contemporary samples
contemp_012<-fread("~/CMH_clumped_outliers_cg_flipped_final.012", na.strings = "-1")
names(contemp_012)<-paste("snp",seq(1:154),sep = "_")
contemp_012<-contemp_012/2 #rescale genotypes to allele frequencies i.e. 0,1,2, to 0,0.5,1
inds<-fread("~/CMH_clumped_outliers_cg_toflip.012.indv",header = F)
contemp_012$sample<-as.character(inds$V1)

#import sample metadata
info<-read.table("~/3waymerged_sampleinfo.txt",sep="\t",header = T) #year and env info
head(info)

#merge files
library(dplyr)
herb_metadata_012<-inner_join(herb_012,info,by="sample")
contemp_metadata_012<-inner_join(contemp_012,info,by="sample")
head(contemp_metadata_012)
contemp_metadata_012$year<-2018
herb_contemp_meta<-rbind(herb_metadata_012,contemp_metadata_012)

#convert from wide to long format
long_012<-melt(herb_contemp_meta,id.vars=c("sample","env","sex","lat","long","year","state"))
long_012 <- long_012 %>% filter(env != "")

#how many samples do we have predating the year 1870?
nrow(herb_metadata_012[herb_metadata_012$year < 1870,])

######### 
#estimate selection by env
#########

#here ag includes dist (since !=Nat)
ag<-herb_contemp_meta %>% filter(env!="") %>% filter(env!="Nat") %>% filter(year>1969)#note herb_metadata_012 vs herb_contemp_meta
nat<-herb_contemp_meta %>% filter(env!="") %>% filter(env=="Nat") %>% filter(year>1969)

#count of historical samples
nrow(nat %>% filter(year < 2018))
nrow(ag %>% filter(year < 2018))

#initialize variables for collecting output of logistic regression over all alleles (n=154)
pval<-list()
zval<-list()
slope<-list()
afmin<-list()
afmax<-list()
slope_error<-list()
newdata <- data.frame(year = c(1870, 2018)) #min and max years to estimate AF at

for (i in 1:154) {
  test<-ag[,..i] #allele number i (since matrix of genos)
  test$year<-ag$year
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial") #regress genotypes on year to get predicted allele frequency and estimate of selection
  test<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test$`Pr(>|z|)`[2] #significance of selection
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2] #selection (haploid currently, multiply by 2 for diploid assuming additivity)
  slope_error[[i]]<-test$`Std. Error`[2] #selection SE
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1] #predicted AF at start of time series
  afmax[[i]]<-probabilities[2] #predicted AF at end of time series
  
}

#merge each vector output into dataframe of results
results_ag<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error),afmin=unlist(afmin),afmax=unlist(afmax))
results_ag$variable<-paste("snp",seq(1:154),sep = "_")

#some summaries
mean(results_ag$afmax) #AF at 2018
mean(results_ag$afmax-results_ag$afmin) #AF change
min(results_ag[results_ag$slope < 1 & results_ag$slope > -1,]$slope)*2 #min selection
max(results_ag[results_ag$slope < 1 & results_ag$slope > -1,]$slope)*2 #max selection


#do the same for individuals in natural environments
year=nat$year
for (i in 1:154) {
  geno<-nat[[i]]
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

#natural results summary
results_nat<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error),afmin=unlist(afmin),afmax=unlist(afmax))
results_nat$variable<-paste("snp",seq(1:154),sep = "_")
pos<-results_nat[results_nat$slope > 0,]
mean(results_nat$afmax-results_nat$afmin)

#working with long format here
long_nat<-long_012 %>% filter(env=="Nat")
long_012_nat<-inner_join(long_nat,results_nat,by="variable")
long_012_nat$pos<-long_012_nat$slope>0
long_012_nat$diff<-long_012_nat$afmax-long_012_nat$afmin
mean(long_012_nat$diff) #mean change
min(results_nat[results_nat$slope < 1 & results_nat$slope > -1,]$slope)*2 #min selection
max(results_nat[results_nat$slope < 1 & results_nat$slope > -1,]$slope)*2 #max selection

long_ag<-long_012 %>% filter(env!="Nat" & env != "")
long_012_ag<-inner_join(long_ag,results_ag,by="variable")
long_012_ag$pos<-long_012_ag$slope>0
long_012_ag$diff<-long_012_ag$afmax-long_012_ag$afmin
mean(long_012_ag$diff)

long_012_both<-rbind(long_012_nat,long_012_ag)

##########
#Joint estimate of s, using logistic regression approach
###########

long_012_both_pregreen<-long_012_both[long_012_both$year < 1960,]
long_012_both_postgreen<-long_012_both[long_012_both$year > 1960,]

long_012_both_postgreen_ag<-long_012_both_postgreen[long_012_both_postgreen$env != "Nat",]
long_012_both_pregreen_ag<-long_012_both_pregreen[long_012_both_pregreen$env != "Nat",]

long_012_both_postgreen_nat<-long_012_both_postgreen[long_012_both_postgreen$env == "Nat",]
long_012_both_pregreen_nat<-long_012_both_pregreen[long_012_both_pregreen$env == "Nat",]

old_timespan <- data.frame(year = c(1870, 1870, 1960, 1960), env=c("Ag","Nat","Ag","Nat") )
new_timespan <- data.frame(year = c(1960, 1960, 2018, 2018), env=c("Ag","Nat","Ag","Nat") )

#selection over the full timespan (1870-2018) in agricultural environments
fulltime_ag<-glm(data=long_ag, value ~ year + variable, family="binomial")
s<-summary(fulltime_ag) 
s$coefficients[[2]] * 2 #since diploid selection

#selection over the full timespan (1870-2018) in natural environments
fulltime_nat<-glm(data=long_nat, value ~ year + variable, family="binomial")
s<-summary(fulltime_nat)
s$coefficients[[2]] * 2 #since diploid selection


#pre and post green revolution

#1870-1960 in agricultural environments
bothenv_pre_ag<-glm(data=long_012_both_pregreen_ag, value ~ year + variable, family="binomial")
s<-summary(bothenv_pre_ag) 
s$coefficients[[2]] * 2 #since diploid selection

#1870-1960 in natural environments
bothenv_pre_nat<-glm(data=long_012_both_pregreen_nat, value ~ year + variable, family="binomial")
#Anova(bothenv_pre_ag,type=3)
s<-summary(bothenv_pre_nat) 
s$coefficients[[2]] * 2 #since diploid selection


#1960-2018 in agricultural environments
bothenv_post_ag<-glm(data=long_012_both_postgreen_ag, value ~ year + variable, family="binomial")
#Anova(bothenv_post_ag,type=3)
s<-summary(bothenv_post_ag)
s$coefficients[[2]] * 2 

#1960-2018 in agricultural environments
bothenv_post_nat<-glm(data=long_012_both_postgreen_nat, value ~ year + variable, family="binomial")
#Anova(bothenv_post_ag,type=3)
s<-summary(bothenv_post_nat)
s$coefficients[[2]] * 2 

#################
#estimate AF at given time point, i.e. 1900 as presented in main text
#################

newdata <- data.frame(year = c(1875,1900,1925,1950,1975,2000,2010,2018))

af1<-list()
af2<-list()
af3<-list()
af4<-list()
af5<-list()
af6<-list()
af7<-list()
af8<-list()



for (i in 1:154) {
  test<-nat[,..i]
  test$year<-nat$year
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  af1[[i]]<-probabilities[1]
  af2[[i]]<-probabilities[2]
  af3[[i]]<-probabilities[3]
  af4[[i]]<-probabilities[4]
  af5[[i]]<-probabilities[5]
  af6[[i]]<-probabilities[6]
  af7[[i]]<-probabilities[7]
  af8[[i]]<-probabilities[8]
  
}

results_nat_moreprob<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),
                                 af1=unlist(af1),af2=unlist(af2),af3=unlist(af3),af4=unlist(af4),
                                 af5=unlist(af5),af6=unlist(af6),af7=unlist(af7),af8=unlist(af8))


ag<-herb_contemp_meta %>% filter(env!="") %>% filter(env!="Nat")
for (i in 1:154) {
  test<-ag[,..i]
  test$year<-ag$year
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  af1[[i]]<-probabilities[1]
  af2[[i]]<-probabilities[2]
  af3[[i]]<-probabilities[3]
  af4[[i]]<-probabilities[4]
  af5[[i]]<-probabilities[5]
  af6[[i]]<-probabilities[6]
  af7[[i]]<-probabilities[7]
  af8[[i]]<-probabilities[8]
  
  
}

results_ag_moreprob<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),
                                af1=unlist(af1),af2=unlist(af2),af3=unlist(af3),af4=unlist(af4),
                                af5=unlist(af5),af6=unlist(af6),af7=unlist(af7),af8=unlist(af8))

std <- function(x) sd(x)/sqrt(length(x))
results_ag_moreprob
ag_means<-melt(colMeans(results_ag_moreprob[,c(4:11)]))
ag_means$year<-c(1875,1900,1925,1950,1975,2000,2010,2018)
ag_means$se<-unlist((sapply(results_ag_moreprob[,c(4:11)],std)))
ag_means$env<-"Ag"

nat_means<-melt(colMeans(results_nat_moreprob[,c(4:11)]))
nat_means$year<-c(1875,1900,1925,1950,1975,2000,2010,2018)
nat_means$se<-unlist(sapply(results_nat_moreprob[,c(4:11)],std))
nat_means$env<-"Nat"
env_means<-rbind(ag_means,nat_means)


###################
#Segmented logistic regression model comparison
###################

library(segmented)
long_012_both_com<-long_012_both[complete.cases(long_012_both),]
out.lm <- glm(value ~ year * env, family="binomial", data = long_012_both_com)
summary(out.lm)

my.seg.k1 <- segmented(out.lm, seg.Z = ~ year, psi=1950,
                       control = seg.control(display = FALSE)
)
summary(my.seg.k1)

my.seg.k2 <- segmented(out.lm, seg.Z = ~ year, psi=list(year=c(1900,2000)),
                       control = seg.control(display = FALSE)
)
summary(my.seg.k2)

AIC(my.seg.k1,my.seg.k2,out.lm)
# get the breakpoints
my.seg.k2$psi #best model is one with two breakpoints, the second falling in 1961



#############
#figure 3b
##############

ggplot() +
  geom_histogram(data=results_nat,aes(slope*2),alpha=.5, color="black",fill="black",bins = 50) +
  geom_histogram(data=results_ag,aes(slope*2),alpha=.6,color="white",fill="grey",bins = 50) +
  theme_classic() +
  xlab("Strength of Selection") +
  ylab("Count") +
  xlim(-0.2,0.4)

#############
#figure S5
#############

results_nat$env<-"Natural"
results_ag$env<-"Agricultural"
median(results_nat$slope)
median(results_ag$slope)

both_sel<-rbind(results_nat[,c("slope","slope_error","variable","env")],results_ag[,c("slope","slope_error","variable","env")])

both_sel2<-both_sel[duplicated(both_sel$variable) | duplicated(both_sel$variable, fromLast=TRUE),]
hist(both_sel2$slope_error,breaks=100)
both_sel2<-both_sel2[both_sel2$slope_error < 1,]
library(tidyr)

both_sel2 %>%
  ggplot() +
  geom_ribbon(aes(ymin=(slope*2)-(slope_error*2),ymax=(slope*2)+(slope_error*2),group=variable, x=env), alpha=.1) +
  geom_point( aes(env,slope*2,group=variable),alpha=.5) +
  geom_line( aes(env,slope*2,group=variable),alpha=.3) +
  theme_bw() +
  labs(y="Strength of Selection",x="Environment") +
  coord_cartesian(ylim=c(-.3,.3))

#############
#plot figure 3C and fig S7
#############


library()
hist_samps<-long_012_both %>% filter(year < 2018) #%>% filter(env=="Ag" | env == "Nat")
contemp_samps<-long_012_both %>% filter(year == 2018) %>% filter(env=="Ag" | env == "Nat")
contemp_samps_thin <-  contemp_samps[-seq(2, NROW(contemp_samps), by = 2),]
contemp_samps_thin <-  contemp_samps_thin[-seq(2, NROW(contemp_samps_thin), by = 2),]
contemp_samps_thin <-  contemp_samps_thin[-seq(2, NROW(contemp_samps_thin), by = 2),]

matched<-rbind(hist_samps,contemp_samps)
matched[matched$env=="Dist"]$env <- "Ag"

avg_n750 <- long_012_both2 %>% filter(year >1869) %>%  #filter(year < 2018) %>%
  #filter(env=="Ag" | env == "Nat") %>%
  arrange(year,variable) %>%
  mutate(id = row_number()) %>% mutate(win=floor(id/1000)) %>% 
  group_by(env, win) %>% summarise(mean_af=mean(value,na.rm=T), sd.mpg = sd(value, na.rm = TRUE),
                                   n.mpg = n(), mean_year=mean(year)) %>%
  mutate(se.mpg = sd.mpg / sqrt(n.mpg),
         lower.ci.mpg = mean_af - qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg,
         upper.ci.mpg = mean_af + qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg)


names(avg_n750)<-c("env","win","value","sd.af","n.af","year","se.af","lower.ci.af","upper.ci.af")
avg_n750_old<-avg_n750[avg_n750$year < 2018,]

avg_n750_new <- long_012_both2 %>% filter(year == 2018) %>%  #filter(year < 2018) %>%
  #filter(env=="Ag" | env == "Nat") %>%
  arrange(year,variable) %>%
  mutate(id = row_number()) %>% mutate(win=floor(id/1000)) %>% 
  group_by(env) %>% summarise(mean_af=mean(value,na.rm=T), sd.mpg = sd(value, na.rm = TRUE),
                              n.mpg = n(), mean_year=mean(year)) %>%
  mutate(se.mpg = sd.mpg / sqrt(n.mpg),
         lower.ci.mpg = mean_af - qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg,
         upper.ci.mpg = mean_af + qt(1 - (0.05 / 2), n.mpg - 1) * se.mpg)

names(avg_n750_new)<-c("env","value","sd.af","n.af","year","se.af","lower.ci.af","upper.ci.af")

mean_by_nobs_year <- 
  long_012_both %>% filter(year >1869) %>%  #filter(year < 2018) %>%
  filter(env=="Ag" | env == "Nat") %>%
  arrange(year,variable) %>%
  mutate(id = row_number()) %>% mutate(win=floor(id/750)) %>% 
  group_by(env, win) %>% dplyr:::summarise(value=mean(value,na.rm=T), year=mean(year))

library(patchwork)
anth<-read.csv("~/Downloads/total-agricultural-land-use-per-person.csv")
head(anth)
library(dplyr)
anth_NA<-anth %>% filter(Entity== "Canada" | Entity == "United States")  %>% group_by(Year) %>% summarise(cropuse=sum(Agricultural.land.per.capita..HYDE..2017..)) 
dens1 <- 
  anth_NA %>% filter(Year > 1860) %>%
  ggplot(aes(x = Year,y=cropuse-2.96026339)) + #standardize by contemp  
  geom_area(color = "black", fill="grey60",stat="identity") + 
  theme_void() +
  scale_color_discrete() +
  theme(
    plot.margin = margin(0, 0,0,0, "cm"))


fig3A<-matched %>%
  #filter(year<2018) %>%
  ggplot(aes(year,value, color=env)) +
  geom_jitter(height = .05,alpha=.10,cex=.75) +
  #geom_smooth(method = "loess", 
  #            method.args = list(span = 0.7), 
  #            se = T, data=test[test$year < 1913,]) +
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), se = T, data=matched[matched$year < 1960 & matched$year > 1870,]) +
  # geom_line(stat="smooth",method = "glm", 
  #           method.args = list(family = "binomial"), se = T, data=matched[matched$year < 1912 & matched$year > 1870,]) +
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), se = T, data=matched[matched$year < 2020 & matched$year > 1960,]) +
  
  theme_bw() +
  geom_point(data=avg_n750_old,cex=3,alpha=.8) +
  geom_errorbar(data=avg_n750_old, aes(ymax=upper.ci.af, ymin=lower.ci.af, x=year),width=0) +
  geom_point(data=avg_n750_new,cex=6,alpha=.8) +
  geom_errorbar(data=avg_n750_new, aes(ymax=upper.ci.af, ymin=lower.ci.af, x=year),width=0) +  #geom_errorbar(data=mean_by_nobs_year, aes(ymin=lower.ci.mpg,ymax=upper.ci.mpg, x=year, color=env),width=2) +
  scale_color_manual(values=c("grey","grey40")) +
  labs(y="Allele Frequency",x="Year") +
  geom_vline(xintercept = 1960,lty="dashed") +
  # xlim(1870,2020) +
  geom_rug(data=matched, aes(x = year,fill=env), outside = T,color="black") +
  theme(axis.text.x = element_text(vjust=-2),
        axis.title.x = element_text(vjust=-1.5),
        axis.ticks.length.x = unit(.4, "cm"),
        axis.ticks.length.y = unit(.2, "cm"),
        axis.ticks.x = element_line(color="lightgrey")) +
  coord_cartesian(clip="off",xlim=c(1865,2020)) 

dens1 + plot_spacer() + fig3A +
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(6,0),
    heights = c(1, 4)
  ) 


fig_S7<-long_012_both %>% filter(year >1869) %>% filter(env=="Ag" | env == "Nat") %>%
  ggplot(aes(year,value, color=env)) + 
  geom_jitter(height = .05,alpha=.1,cex=1) +
  stat_smooth(method = "lm", formula = y ~ splines::bs(x, 3),level=.95) +
  theme_bw() +
  geom_point(data=avg_n750_old,cex=3,alpha=.8) +
  geom_errorbar(data=avg_n750_old, aes(ymax=upper.ci.af, ymin=lower.ci.af, x=year),width=0) +
  geom_point(data=avg_n750_new,cex=6,alpha=.8) +
  geom_errorbar(data=avg_n750_new, aes(ymax=upper.ci.af, ymin=lower.ci.af, x=year),width=0) +
  scale_color_manual(values=c("grey","grey40")) +
  labs(y="Ag-Allele Frequency",x="Year") 


#########
#Figure S4 and contemporary frequencies by env
########
freqs<-contemp_metadata_012 %>% filter(env!="Dist") %>% group_by(env) %>% 
  summarise_at(.vars = paste("snp",1:154,sep="_"),.funs = list(mean=mean), na.rm=T)

freqs_t<-as.data.frame(t(freqs))
names(freqs_t)<-c("Ag","Nat")
freqs_t<-freqs_t[-1,]
freqs_t$Ag<-as.numeric(as.character(freqs_t$Ag))
freqs_t$Nat<-as.numeric(as.character(freqs_t$Nat))

p<-ggplot(data=freqs_t, aes(Ag, Nat)) +
  #geom_rug(col=rgb(.5,0,0,alpha=.8),) +
  annotate("rect",ymin = min(freqs_t$Nat), ymax=max(freqs_t$Nat), xmin=0, xmax=1,alpha=.15) +
  annotate("rect",xmin = min(freqs_t$Ag), xmax=max(freqs_t$Ag), ymin=0, ymax=1,alpha=.15) +
  geom_point() +
  ylim(0,1) +
  xlim(0,1) +
  theme_bw() +
  labs(y="Nat Frequency",x="Ag Frequency")

ggExtra::ggMarginal(p, type = "histogram",size=3)


mean(freqs_t$Ag)
mean(freqs_t$Nat)

##########
#figure 3A
##########


long_012_ag
long_012_nat
quantile(results_ag$diff)

long_012_nat$quant<-"<89%"
long_012_nat[long_012_nat$diff <= quantile(results_ag$diff)[2],]$quant <-"<2%"
long_012_nat[long_012_nat$diff <= quantile(results_ag$diff)[3] & long_012_nat$diff > quantile(results_ag$diff)[2],]$quant <-"<19%"
long_012_nat[long_012_nat$diff <= quantile(results_ag$diff)[4] & long_012_nat$diff > quantile(results_ag$diff)[3],]$quant <-"<48%"
long_012_nat$quant<-as.factor(long_012_nat$quant)
long_012_nat$quant <- factor(long_012_nat$quant, levels = c("<2%","<19%", "<48%", "<89%"))


long_012_ag$quant<-"<89%"
long_012_ag[long_012_ag$diff <= quantile(results_ag$diff)[2],]$quant <-"<2%"
long_012_ag[long_012_ag$diff <= quantile(results_ag$diff)[3] & long_012_ag$diff > quantile(results_ag$diff)[2],]$quant <-"<19%"
long_012_ag[long_012_ag$diff <= quantile(results_ag$diff)[4] & long_012_ag$diff > quantile(results_ag$diff)[3],]$quant <-"<48%"
long_012_ag$quant<-as.factor(long_012_ag$quant)
long_012_ag$quant <- factor(long_012_ag$quant, levels = c("<2%","<19%", "<48%", "<89%"))


long_012_nat$sig<-"yes"
long_012_nat[long_012_nat$pval > 0.05,]$sig<-"no"

nat_plot <- long_012_nat %>% filter(year>1870) %>%
  ggplot(aes(year,value,group=variable,color=quant,alpha=sig)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F,cex=1.5) +
  theme_bw() +
  labs(y="Allele Frequency",x="Year",color="Quantile") +
  scale_color_viridis_d(direction = -1,end = .9, option = "C") +
  #scale_linetype_manual(values = c("21","solid")) +
  scale_alpha_manual(values =c(0.3, 1)) +
  geom_rug(aes(y=NULL),outside = T,color="black") +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(vjust=-1),
        axis.ticks.length.x =unit(.4, "cm"),
        axis.ticks.x = element_line(colour="lightgrey"))

long_012_ag$sig<-"yes"
long_012_ag[long_012_ag$pval > 0.05,]$sig<-"no"

ag_plot<-long_012_ag %>% filter(year>1870) %>%
  ggplot(aes(year,value,group=variable,color=quant,alpha=sig)) + 
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F,cex=1.5) +
  theme_bw() +
  labs(y="Allele Frequency",x="Year",color="Quantile") +
  scale_color_viridis_d(direction =-1,end = .9, option = "C") +
  theme(legend.position = "none") +
  scale_linetype_manual(values = c("21","solid")) +
  scale_alpha_manual(values =c(0.3,1)) +
  geom_rug(aes(y=NULL),outside = T,color="black") +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(vjust=-1),
        axis.ticks.length.x =unit(.4, "cm"),
        axis.ticks.x = element_line(colour="lightgrey"))

library(cowplot)
library(gridExtra)

grid.arrange(ag_plot,nat_plot,ncol=2,widths=c(0.46,0.54))


#############
#figure S8
#############

ag_old <- ag %>% filter(year < 1960)
nat_old<- nat %>% filter (year < 1960)
ag_recent <- ag %>% filter(year > 1960)
nat_recent <- nat %>% filter(year > 1960)

slope_error<-list()
for (i in 1:154) {
  test<-ag_old[,..i]
  test$year<-ag_old$year
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  slope_error[[i]]<-test$`Std. Error`[2]
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1]
  afmax[[i]]<-probabilities[2]
  
}

results_ag_old<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error))
results_ag_old$variable<-paste("snp",seq(1:154),sep = "_")
results_ag_old <- results_ag_old %>% filter(slope > -1 & slope < 1)
results_ag_old <- results_ag_old %>% filter(slope_error < 1)

mean(results_ag_old$slope)
mean(results_ag_old$slope_error)


for (i in 1:154) {
  test<-ag_recent[,..i]
  test$year<-ag_recent$year
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  slope_error[[i]]<-test$`Std. Error`[2]
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1]
  afmax[[i]]<-probabilities[2]
  
}
results_ag_recent<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error))
results_ag_recent$variable<-paste("snp",seq(1:154),sep = "_")
results_ag_recent <- results_ag_recent %>% filter(slope > -1 & slope < 1)
results_ag_recent <- results_ag_recent %>% filter(slope_error < 1)

mean(results_ag_recent$slope)
mean(results_ag_recent$slope_error)


slope_error<-list()
for (i in 1:154) {
  test<-nat_old[,..i]
  test$year<-nat_old$year
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  slope_error[[i]]<-test$`Std. Error`[2]
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1]
  afmax[[i]]<-probabilities[2]
  
}
results_nat_old<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error))
results_nat_old$variable<-paste("snp",seq(1:154),sep = "_")
results_nat_old <- results_nat_old %>% filter(slope > -1 & slope < 1)
results_nat_old <- results_nat_old %>% filter(slope_error < 1)

mean(results_nat_old$slope)
mean(results_nat_old$slope_error)


for (i in 1:154) {
  test<-nat_recent[,..i]
  test$year<-nat_recent$year
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial")
  test<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test$`Pr(>|z|)`[2]
  zval[[i]]<-test$`z value`[2]
  slope[[i]]<-test$Estimate[2]
  slope_error[[i]]<-test$`Std. Error`[2]
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1]
  afmax[[i]]<-probabilities[2]
  
}
results_nat_recent<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),slope_error=unlist(slope_error))
results_nat_recent$variable<-paste("snp",seq(1:154),sep = "_")
results_nat_recent <- results_nat_recent %>% filter(slope > -1 & slope < 1)
results_nat_recent <- results_nat_recent %>% filter(slope_error < 1)

mean(results_nat_recent$slope)
mean(results_nat_recent$slope_error)


old_s_hist<-ggplot() +
  geom_histogram(data=results_nat_old,aes(slope*2),alpha=.5, color="black",fill="black",bins = 50) +
  geom_histogram(data=results_ag_old,aes(slope*2),alpha=.6,color="white",fill="grey",bins = 50) +
  theme_classic() +
  xlab("Strength of Selection") +
  ylab("Count") +
  xlim(-1,1)

recent_s_hist<-ggplot() +
  geom_histogram(data=results_nat_recent,aes(slope*2),alpha=.5, color="black",fill="black",bins = 50) +
  geom_histogram(data=results_ag_recent,aes(slope*2),alpha=.6,color="white",fill="grey",bins = 50) +
  theme_classic() +
  xlab("Strength of Selection") +
  ylab("Count") +
  xlim(-1,1)

grid.arrange(old_s_hist, recent_s_hist, ncol=2)


#####################
#figure S6
######################


wrange<-read.table("cmh_maf01_clumped_wranges.gff.fdr.final.match",sep="\t")
head(wrange)
names(wrange)<-c("scaf","x","pos","ass_pval","n","n_sig","length","scafz","xx","y","z","zz","zzz","zzzz","zzzzz","ann")
head(wrange)


wrange$variable<-paste("snp_",seq(1:154),sep="")

ag_wrange<-inner_join(wrange,results_ag,by="variable")
names(ag_wrange)

ag_wrange$totalchange<-ag_wrange$afmax-ag_wrange$afmin
ag_range_wxpehh_wrecomb<-read.table("xpehh_aglogistic_recomb.bed",sep="\t")
#head(ag_range_wxpehh_wrecomb)
names(ag_range_wxpehh_wrecomb)[c(1,3,4,5,6,7,17,18,19,21,30,32,36,40)]<-c("scaf","pos","p_val","total_n","sig_n","length", "variable",
                                                                          "change_pval","slope","start_frq","end_fq","total_change","xpehh","recomb")
ag_wrange_all<-inner_join(ag_wrange,ag_range_wxpehh_wrecomb[,c(17,36,40)])
hist(ag_wrange_all$slope)
ag_wrange_all<-ag_wrange_all[ag_wrange_all$slope_error < 50,]

lm1<-lm(data=ag_wrange_all, slope ~ n + length + ass_pval + xpehh + recomb )
summary(lm1)
Anova(lm1, type=3)

lm2<-lm(data=ag_wrange_all, totalchange ~ n + length + ass_pval + xpehh + recomb )
summary(lm2)
Anova(lm2, type=3)

library(lsmeans)
lsmean_slope<- lsmip(lm1, ~ n, ylab = "Total Change",
                     at=list(n=c(0, 100, 200, 300, 400, 500, 600,700)),
                     xlab="# SNPs in clump (r2 > 0.25)",
                     type="response", plot=F)



lsmean_totalchange<- lsmip(lm2, ~ n, ylab = "Total Change",
                    at=list(n=c(0, 100, 200, 300, 400, 500, 600,700)),
                    xlab="# SNPs in clump (r2 > 0.25)",
                    type="response", plot=F)


ggplot() +
  geom_point(data=ag_wrange_all, aes(n, slope*2 )) +
  geom_smooth(data=lsmean_slope, aes(n, yvar*2), method="lm",color="black") +
  geom_ribbon(data=lsmean_slope, aes(ymax=yvar*2+SE*2,ymin=yvar*2-SE*2,x=n),alpha=.2) +
  theme_bw() +
  labs(y="Selection", x="SNPs in linkage with focal allele \nat r2 > 0.25") 
  
ggplot() +
    geom_point(data=ag_wrange_all, aes(n, totalchange )) +
    geom_smooth(data=lsmean_totalchange, aes(n, yvar), method="lm",color="black") +
    geom_ribbon(data=lsmean_totalchange, aes(ymax=yvar+SE,ymin=yvar-SE,x=n),alpha=.2) +
    theme_bw() +
    labs(y="Total AF Change", x="SNPs in linkage with focal allele \nat r2 > 0.25") +
    coord_cartesian(ylim = c(-.5, 1)) 

