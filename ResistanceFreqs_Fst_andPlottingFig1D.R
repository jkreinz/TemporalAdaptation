library(ggplot2)
library(dplyr)

############
#use 012 coding in the paper
commongarden<-read.csv("~/Downloads/Common Garden Seq Stats - moremeta (1).csv")
commongarden$newcoding<-0
commongarden$newcoding[commongarden$EPSPS_scaledcopy > 1.5 & commongarden$EPSPS_scaledcopy < 2.5] <- 1
commongarden$newcoding[commongarden$EPSPS_scaledcopy > 2.5] <- 2
commongarden$newcoding

#differentiation for EPSPS amplification?
commongarden %>% group_by(Env.x) %>% dplyr:::summarise(mean_EPSPScopy=mean(EPSPS_scaledcopy))
pop_mean<-commongarden %>% group_by(Pair,Env.x) %>% dplyr:::summarise(mean_EPSPScopy=mean(EPSPS_scaledcopy))
Anova(lm(data=commongarden, EPSPS_scaledcopy ~ Env.x + Pair),type = 3)

#calculate Fst
p1s<-list()
p2s<-list()
n1s<-list()
n2s<-list()
fsts<-list()

for (i in sort(unique(commongarden$Pair))) {
  pair<-commongarden[commongarden$Pair == i,]
  ag<-pair[pair$Env.x=="Ag",]
  nat<-pair[pair$Env.x=="Nat",]
  
  p1<-mean((ag$newcoding))/2
  p2<-mean((nat$newcoding))/2
  p1s[i]<-p1
  p2s[i]<-p2
  
  n1<-nrow(ag)*2
  n2<-nrow(nat)*2
  n1s[i]<-n1
  n2s[i]<-n2
  
  num=((p1-p2)^2) - ((p1*(1-p1))/(n1-1)) - ((p2*(1-p2))/(n2-1))
  denom=p1*(1-p2) + p2*(1-p1)
  
  fsts[i]<-num/denom
}

unlist(p1s)
unlist(p2s)

mean(unlist(p1s),na.rm=T)
mean(unlist(p2s),na.rm=T)

fst_bypair<-data.frame(ag_p=unlist(p1s),nat_p=unlist(p2s),ag_n=unlist(n1s), nat_n=unlist(n2s),fst=unlist(fsts))
(mean(fst_bypair$fst,na.rm = T)) #within pair mean Fst for EPSPS amplification


#now across envs, regardless of pair

commongarden$newcoding<-0
commongarden$newcoding[commongarden$EPSPS_scaledcopy > 1.5 & commongarden$EPSPS_scaledcopy < 2.5] <- 1
commongarden$newcoding[commongarden$EPSPS_scaledcopy > 2.5] <- 2
commongarden$newcoding

ag<-commongarden[commongarden$Env == "Ag",]
nat<-commongarden[commongarden$Env == "Nat",]

(p1<-mean((ag$newcoding))/2)
(p2<-mean((nat$newcoding))/2)
n1<-nrow(ag)*2
n2<-nrow(nat)*2

num=((p1-p2)^2) - ((p1*(1-p1))/(n1-1)) - ((p2*(1-p2))/(n2-1))
denom=p1*(1-p2) + p2*(1-p1)

((hudsons_fst=num/denom)) #across env Fst


#############
#for all other resistance alleles
#############


tsr<-read.table("~/TSRgenos_commongarden.txt",header=T)
head(tsr)
#meta<-read.table("herbcontemp_poplist.txt",na.strings = c("","NA","?"),sep="\t",header = F)
meta<-read.table("~/contemp_wpair.txt",na.strings=".")
names(meta)<-c("Sample","Env","Sex","Pair")
library(car)
library(dplyr)
tsr_wmeta<-inner_join(tsr,meta, by="Sample")

#check for sig differentiation across environments
Anova(lm(data=tsr_wmeta, ALS653 ~ Env + Pair),type = 3)
Anova(lm(data=tsr_wmeta, ALS574 ~ Env + Pair),type = 3)
Anova(lm(data=tsr_wmeta, ALS197 ~ Env + Pair),type = 3)
Anova(lm(data=tsr_wmeta, ALS122 ~ Env + Pair),type = 3)
Anova(lm(data=tsr_wmeta, ALS376 ~ Env + Pair),type = 3)
Anova(lm(data=tsr_wmeta, EPSPS210 ~ Env + Pair),type = 3)
Anova(lm(data=tsr_wmeta, PPOdel ~ Env + Pair),type = 3)

#get frequencies
ag<-tsr_wmeta[tsr_wmeta$Env == "Ag",]
nat<-tsr_wmeta[tsr_wmeta$Env == "Nat",]

coords<-read.table("contemp_coords.txt")
names(coords)<-c("Sample","long","lat")
tsr_wmeta<-inner_join(tsr_wmeta, coords, by="Sample")

bypair<-tsr_wmeta %>% group_by(Pair, Env) %>% 
  group_by(Pair, Env, N = n(), add = TRUE) %>% 
  summarise(n = n(), ALS653=mean(ALS653)/2 , ALS197=mean(ALS197)/2, ALS122=mean(ALS122)/2, 
            ALS376 = mean(ALS376)/2, ALS574=mean(ALS574)/2, EPSPS210=mean(EPSPS210)/2, 
            PPOdel=mean(PPOdel)/2, PPOdel_std=std(na.omit(PPOdel)))

byenv_bypair<-bypair %>% group_by(Env) %>% 
  dplyr:::summarise(n = n(), ALS653=mean(ALS653), ALS197=mean(ALS197), ALS122=mean(ALS122), 
                    ALS376 = mean(ALS376), ALS574=mean(ALS574), EPSPS210=mean(EPSPS210), 
                    PPOdel=mean(PPOdel))

byenv_bypair<-as.data.frame(t(byenv_bypair[-c(1:2)])) #freqs of TSR in nat vs ag
names(byenv_bypair)<-c("ag","nat")
byenv_bypair$ag-byenv_bypair$nat #diff in freq by env


####
#fst for PPOdel
####

(p1<-mean((ag$PPOdel))/2)
(p2<-mean((nat$PPOdel))/2)
n1<-nrow(ag)*2
n2<-nrow(nat)*2

num=((p1-p2)^2) - ((p1*(1-p1))/(n1-1)) - ((p2*(1-p2))/(n2-1))
denom=p1*(1-p2) + p2*(1-p1)

((hudsons_fst=num/denom)) #across env Fst for PPO

################
#plot Figure 2D
################

std <- function(x) sd(x)/sqrt(length(x))
bypair<-tsr_wmeta %>% group_by(Pair, Env) %>% 
  group_by(Pair, Env, N = n(), add = TRUE) %>% 
  summarise(n = n(), ALS653=mean(ALS653)/2 , ALS197=mean(ALS197)/2, ALS122=mean(ALS122)/2, 
            ALS376 = mean(ALS376)/2, ALS574=mean(ALS574)/2, EPSPS210=mean(EPSPS210)/2, 
            PPOdel=mean(PPOdel)/2, PPOdel_std=std(na.omit(PPOdel)))

byenv<-tsr_wmeta %>% 
  group_by(Env, N = n(), add = TRUE) %>% 
  summarise(n = n(), ALS653=mean(ALS653)/2 , ALS197=mean(ALS197)/2, ALS122=mean(ALS122)/2, 
            ALS376 = mean(ALS376)/2, ALS574=mean(ALS574)/2, EPSPS210=mean(EPSPS210)/2, 
            PPOdel=mean(PPOdel)/2, PPOdel_std=std(na.omit(PPOdel)))


byenv_long<-as.data.frame(t(byenv[,4:10]))
names(byenv_long)<-c("Ag","Nat")
mean(byenv_long$Ag - byenv_long$Nat)

ymin_e<-bypair$PPOdel-bypair$PPOdel_std
ymax_e<-bypair$PPOdel+bypair$PPOdel_std


PPO<-ggplot(bypair, aes(Env,PPOdel, group=Pair)) +
  geom_point() +
  geom_line() +
  theme_bw() +  ylim(0,1) +
  ylab("Frequency of PPOdel") 

A574<-ggplot(bypair, aes(Env,ALS574, group=Pair)) +
  geom_point() +
  geom_line() +
  theme_bw() +  ylim(0,1) +
  ylab("Frequency of ALS574")

A653<-ggplot(bypair, aes(Env,ALS653, group=Pair)) +
  geom_point() +
  geom_line() +
  theme_bw() +  ylim(0,1) +
  ylab("Frequency of ALS653")

A376<-ggplot(bypair, aes(Env,ALS376, group=Pair)) +
  geom_point() +
  geom_line() +
  theme_bw() +   ylim(0,1) +
  ylab("Frequency of ALS376")

fst_bypair<-data.frame(ag_p=unlist(p1s),nat_p=unlist(p2s))
names(fst_bypair)<-c("Ag","Nat")
fst_bypair$pair<-seq(1:17)
fst_long<-melt(fst_bypair,id.vars="pair")
EP210<-ggplot(bypair, aes(Env,EPSPS210, group=Pair)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylim(0,1) +
  ylab("Frequency of EPSPS106")

head(fst_bypair)
names(fst_bypair)[1:2]<-c("Ag","Nat")
fst_bypair$pair<-seq(1:17)
fst_long<-melt(fst_bypair[,c(1:2,6)],id.vars="pair")

#check for sig diff at EPSPSamp
Anova(lm(data=fst_long, value ~ variable + pair),type = 3)

EPSPSamp<-ggplot(fst_long, aes(variable,value, group=pair)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylim(0,1) +
  labs(y="Frequency of EPSPS amp.", x="Env")

#finally, figure 2D
grid.arrange(A653, A376, A574, EP210, EPSPSamp, PPO, ncol=6)



