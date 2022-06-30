library(data.table)
herb_012<-fread("~/herb_resSNPs_sort.012",na.strings = "-1")
herb_012<-herb_012[,-1]
#names(herb_012)<-paste("snp",seq(1:7),sep = "_")
names(herb_012)<-c("PPO210","ALS653","ALS574","ALS376","ALS197","ALS122","EPSPS210")
herb_012<-herb_012/2
inds<-fread("~/herb_resSNPs_sort.012.indv",header = F)
herb_012$sample<-gsub("HBO","HB0",inds$V1)



contemp_012<-fread("~/contemporary_resSNPs_sort.012", na.strings = "-1")
contemp_012<-contemp_012[,-1]
names(contemp_012)<-c("PPO210","ALS653","ALS574","ALS376","ALS197","ALS122","EPSPS210")
contemp_012<-contemp_012/2
inds<-fread("~/contemporary_resSNPs_sort.012.indv",header = F)
contemp_012$sample<-as.character(inds$V1)


info<-read.table("~/3waymerged_sampleinfo.txt",sep="\t",header = T)
head(info)


######

herb_metadata_012<-inner_join(herb_012,info,by="sample")
contemp_metadata_012<-inner_join(contemp_012,info,by="sample")
head(contemp_metadata_012)
contemp_metadata_012$year<-2018

herb_metadata_012$ALS574[which(herb_metadata_012$ALS574 == 0.5)] <- 0
herb_metadata_012$ALS574[which(herb_metadata_012$ALS574 == 0.5)]

herb_metadata_012[which(herb_metadata_012$ALS122 == 0.5),]
herb_metadata_012$ALS122[which(herb_metadata_012$ALS122 == 0.5 & herb_metadata_012$year < 2009)]<-0

herb_contemp_meta<-rbind(herb_metadata_012,contemp_metadata_012)

long_012<-melt(herb_contemp_meta,id.vars=c("sample","env","sex","lat","long","year","state"))
long_012 <- long_012 %>% filter(env != "")
#long_012<-melt(herb_metadata_012,id.vars=c("sample","env","sex","lat","long","year","state"))

nrow(herb_metadata_012[herb_metadata_012$year < 1870,])



pval<-list()
zval<-list()
slope<-list()
intcp<-list()
afmin<-list()
afmax<-list()
se<-list()
########
newdata <- data.frame(year = c(1960, 2018))

for (i in 1:7) {
  test<-herb_contemp_meta[,..i]
  test$year<-herb_contemp_meta$year
  test<-test[test$year > 1960, ]
  
  lm1<-glm(data=test, unlist(test[,1]) ~ year, family="binomial")
  test2<-as.data.frame(summary(lm1)$coefficients)
  
  pval[[i]]<-test2$`Pr(>|z|)`[2]
  zval[[i]]<-test2$`z value`[2]
  slope[[i]]<-test2$Estimate[2]
  intcp[[i]]<-test2$Estimate[1]
  se[[i]]<-test2$`Std. Error`[2]
  
  probabilities <- lm1 %>% predict(newdata, type = "response")
  afmin[[i]]<-probabilities[1]
  afmax[[i]]<-probabilities[2]
  
}

results<-data.frame(pval=unlist(pval),zval=unlist(zval),slope=unlist(slope),se=unlist(se),inter=unlist(intcp),afmin=unlist(afmin),afmax=unlist(afmax))
row.names(results)<-c("PPO210","ALS653","ALS574","ALS376","ALS197","ALS122","EPSPS210")
write.table(results,"individualmods_sel_resistancealleles.txt",sep="\t",col.names = T,row.names = T,quote=F)

hist(results$afmax-results$afmin)
mean(results$afmax-results$afmin)
median(results$afmax-results$afmin)
pos<-results[results$slope > 0,]
nrow(pos)
median(pos$afmax-pos$afmin)
neg<-results[results$slope < 0,]
nrow(neg)
median(neg$afmax-neg$afmin)

head(long_012)
tail(long_012)

results$variable<-c("PPO210","ALS653","ALS574","ALS376","ALS197","ALS122","EPSPS210")
long_012_merged<-inner_join(long_012,results,by="variable")
#write.table(long_012,"seven_herbicideresistant_loci_genotypes_throughtime.csv",sep=",",col.names = T,row.names = F)
long_012_merged$pos<-long_012_merged$slope>0
#install.packages("MetBrewer")
library(MetBrewer)
library(ggplot2)


anth<-read.csv("~/Downloads/total-agricultural-land-use-per-person.csv")
#effic<-read.csv("~/Downloads/total-agricultura")

head(anth)

library(dplyr)
anth_NA<-anth %>% filter(Entity== "Canada" | Entity == "United States")  %>% group_by(Year) %>% summarise(cropuse=sum(Agricultural.land.per.capita..HYDE..2017..))


dens1 <- 
anth_NA %>% filter(Year > 1860) %>%
ggplot(aes(x = Year,y=cropuse-2.69030750)) + 
  geom_area( fill="grey75",stat="identity") + 
  #theme(axis.line=element_blank(),
  #     axis.text.x=element_blank(),
  #     axis.text.y=element_blank(),
  #     axis.title.x=element_blank(),
  #     axis.title.y=element_blank(),
  #     legend.position="none",
  #     panel.background=element_blank(),
  #     panel.border=element_blank(),
  #     panel.grid.major=element_blank(),
  #     panel.grid.minor=element_blank(),
  #     plot.background=element_blank()) +
  #theme(axis.line.x = element_line(size = 2)) +
theme_void()

library(patchwork)

herb_s<-ggplot(data=long_012_merged,aes(year,value,group=variable,color=variable)) + 
  geom_jitter(height = .05,width=0,cex=2) +
  stat_smooth(method = "glm",
            method.args = list(family = "binomial"), 
            se = F, alpha=.1,cex=1.25) +
  theme_bw() +
  labs(y="Allele Frequency",x="Year",color="Resistance\nAllele") +
  scale_color_manual(values=met.brewer("Tiepolo", 8)) +
  xlim(1870,2020) +
  geom_vline(xintercept=1960,lty="dashed") +
  geom_rug(aes(y=NULL),outside = T,color="black") +
  coord_cartesian(clip="off") +
  theme(axis.text.x = element_text(vjust=-1),
        axis.ticks.length.x =unit(.4, "cm"),
        axis.ticks.x = element_line(colour="lightgrey"),
        plot.margin = margin(0, 0,0,0, "cm"))
  

library(patchwork)

##########
#figure 3D
dens1 + plot_spacer() + herb_s +
  plot_layout(
    ncol = 2, 
    nrow = 2, 
    widths = c(6,0),
    heights = c(1, 4)
  ) 

#figure 3D inset
ggplot(data=long_012_merged,aes(year,value,group=variable,color=variable)) + 
  #geom_jitter(height = .1) +
  geom_line(stat="smooth",method = "glm", 
            method.args = list(family = "binomial"), 
            se = F, alpha=.9,cex=1.5,level = 0.2) +
  theme_bw() +
  labs(y="Ag-Allele Frequency",x="Year",color="Resistance\nAllele") +
  scale_color_manual(values=met.brewer("Tiepolo", 8)) +
  coord_cartesian(xlim=c(2000,2020)) 


avg_resistance<-glm(data=long_012[long_012$year > 1960 & long_012$env != "Dist",], value ~ year * env + variable, family="binomial")
s<-summary(avg_resistance)
s$coefficients[[2]] * 2 #estimate of selection across environments since 1960
s$coefficients[[2,4]] #p-value

avg_ag_resistance<-glm(data=long_012[long_012$year > 1960 & long_012$env == "Ag",], value ~ year + variable, family="binomial")
summary(avg_ag_resistance)
s<-summary(avg_ag_resistance)
s$coefficients[[2]] * 2 #estimate of selection in ag since 1960
s$coefficients[[2,4]] #p-value

nat_long<-long_012[long_012$year > 1960 & long_012$env == "Nat",]
avg_nat_resistance<-glm(data=nat_long, value ~ year + variable, family="binomial")
summary(avg_nat_resistance)
s<-summary(avg_nat_resistance)
s$coefficients[[2]] * 2 #estimate of selection in nat since 1960
s$coefficients[[2,4]] #p-value

dist_long<-long_012[long_012$year > 1960 & long_012$env == "Dist",]
avg_dist_resistance<-glm(data=dist_long, value ~ year + variable, family="binomial")
summary(avg_dist_resistance)
s<-summary(avg_dist_resistance)
s$coefficients[[2]] * 2 #estimate of selection in dist since 1960
s$coefficients[[2,4]] #p-value

################
#Figure S2
################

costben<-read.csv("~/Desktop/science_submission/costben CMH outliers.csv",header=T)
head(costben)

somecostben<-data.frame(costben$ag_benefit,costben$nat_cost,costben$snp)
names(somecostben)<-c("ag_benefit","nat_cost","allele")
herb<-data.frame(ag_benefit=c(2.58582, 2.08013,1.30689),nat_cost=c(12.9997,2.88402,1.88974),
                 allele=c("PPOdel","EPSPSamp","ALS574"))

library(MASS) # to access Animals data sets
library(scales) # to access break formatting functions

scaleFUN <- function(x) sprintf("%.2f", x)

ggplot() +
  geom_point(data=somecostben,aes(nat_cost,ag_benefit)) +
  geom_point(data=herb,aes(nat_cost,ag_benefit, color=allele)) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=scaleFUN) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=scaleFUN) +
  theme_bw() +
  geom_abline() +
  labs(y=expression(paste("Agricultural  ", italic(frac("s","m")))), x=expression(paste("Natural  ", italic(frac("s","m")))))
  

mean(costben$cost_ben_ratio,na.rm=T)
median(costben$cost_ben_ratio,na.rm=T)



