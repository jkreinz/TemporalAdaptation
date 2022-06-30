#plotting population structure through time

library(pophelper)
library(dplyr)
library(data.table)

setwd("~/3waymerged_allsamps/") #directory with faststructure results for K2

sfiles <- list.files(pattern = "*2.meanQ")
slist <- sortQ(readQ(sfiles))
slist_1 <- alignK(slist)
slist_final<- mergeQ(slist_1)

# collect aligned files. Default prefix is pop
sfiles <- list.files(pattern = "*3.meanQ")
slist <- sortQ(readQ(sfiles))
slist_1 <- alignK(slist)
slist_final[2] <- mergeQ(slist_1)
slist_final<-as.qlist(slist_final)


#read in metadata

pop<-read.table("~/herbcontemp_poplist.txt",na.strings = c("","NA"),sep="\t")
names(pop)<-c("sample","env","sex","lat","long","year")
pop2<-read.table("~/3waymerged_sampleinfo.txt",sep="\t",header = T)
both_pop<-inner_join(pop2, pop, by="sample") #to fill in missing metadata from first file
head(both_pop)

#make sure order of metadata matches output of faststructure (which used plink files as input, so taking the fam)
ch_order<-read.table("~/3waymerged_allsamps_geno20p_maf05.fam") 
names(ch_order)[1]<-c("sample")
both_pop<-inner_join(ch_order[1], pop2, by="sample")

#add K values to metadata 
both_pop$K1<- slist_final[[1]][["Cluster1"]] #var. rudis
both_pop$K2<- slist_final[[1]][["Cluster2"]] #var. tuberculatus

#rename cols for plotting
names(both_pop)[1:7]<-c("sample","env","sex","lat","long","year","state")
both_pop$state[both_pop$state == "Nat"]<-"Ontario"
both_pop$state[both_pop$state == "Walpole"]<-"Ontario"
both_pop$state[both_pop$state == "Essex"]<-"Ontario"

#add timespans
both_pop$time_span<-"2020"
both_pop$time_span[both_pop$year < 1920]<-"1920"
both_pop$time_span[both_pop$year < 1980 & both_pop$year > 1919]<-"1980"
std <- function(x) sd(x)/sqrt(length(x))


#plotting
detach(package:plyr)
both_pop$time_span
std <- function(x) sd(x)/sqrt(length(x))

both_pop$state <- factor(both_pop$state, levels = c("Ontario","Ohio", "Michigan", "Indiana","Illinois","Missouri","Kansas"))
state_byyear<-both_pop %>% group_by(state,time_span) %>% summarise(meanK1=mean(K1), seK1=std(K1), n=n())
byyear<-both_pop %>% filter(year<2015) %>% group_by(time_span) %>% summarise(meanK1=mean(K1), medianK1=median(K1), seK1=std(K1), n=n())
all<-both_pop %>% filter(year<2015) %>% summarise(meanK1=mean(K1), medianK1=median(K1), seK1=std(K1), n=n())

state_byyear$time_span<-as.numeric(state_byyear$time_span)
both_pop$time_span<-as.factor(both_pop$time_span)

###########
#Figure 4A
###########
ggplot(data=both_pop, aes(long,K1, group=time_span, color=time_span)) + 
  geom_point(alpha=.8) +
  #stat_smooth(method = "lm", formula = y ~ x + I(x^2) + I(x^3), size = 1,alpha=.2) +
  #stat_smooth(method = "gam", size = 1,alpha=.2,span=1) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se =T,alpha=.2, level=.95) +
  theme_bw() +
  labs(color = "Time span", x="Longitude", y=substitute(paste("Proportion var. ",italic("rudis")))) +
  scale_color_viridis_d(end =  .9) +
  coord_cartesian(ylim = c(0,1))

###########
#Figure 4B
###########
library(ggridges)
both_pop$time_span<-as.factor(both_pop$time_span)
both_pop %>% filter(state != "Michigan") %>% filter(state != "Indiana" ) %>% 
  ggplot(aes(x = K1, y = interaction(time_span,state), fill=time_span)) +
  #stat_density_ridges(quantile_lines = TRUE, alpha=.7,show.legend = F,calc_ecdf = T) +
  theme_bw() +
  geom_point(color="grey40",alpha=.5) +
  geom_density_ridges(quantile_lines = T, alpha=.7,from=0, to=1, show.legend = F,calc_ecdf = T) +
  scale_fill_viridis_d(end =  .9) +
  xlab(substitute(paste("Proportion var. ",italic("rudis"), " ancestry"))) +
  labs(y="Time Span by State", fill="Time Span") +
  annotate("rect", xmin = 0, xmax = 0, ymin = -.25, ymax = 17,
           alpha = 0,fill = "white") 
# annotate("rect", xmin = 1, xmax = 1.5, ymin = -.5, ymax = 16,
#         alpha = .8,fill = "white") 

##########
#Figure 4C
##########

both_pop$dataset<-"commongarden"
both_pop[both_pop$year == "2015",]$dataset<-"pnas"
both_pop[both_pop$year < "2015",]$dataset<-"herb"

both_pop[both_pop$env == "Dist",]$env<-"Ag"

justherb<-both_pop %>% filter(dataset == "herb") %>% 
  filter(env != "") %>%
  ggplot(aes(long,K1, color=env)) +
  geom_point(alpha=.7) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = T, alpha=.2) +
  theme_bw() +
  xlim(-96,-74) +
  labs(x="Longitude",y="Proportion var. rudis ancestry",title = "Herbarium") +
  scale_color_manual(values=c("grey","grey40"))

together<-both_pop %>% filter(dataset != "PNAS") %>% 
  filter(env != "") %>%
  ggplot(aes(long,K1, color=env)) +
  geom_point(alpha=.7) +
  #geom_line(data=data.frame(spline(d, n=n*10)))
  
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = T,alpha=.2) +
  theme_bw() +
  xlim(-96,-74) +
  labs(x="Longitude",y="Proportion var. rudis ancestry", title = "Herbarium + Contemporary") +
  scale_color_manual(values=c("grey","grey40"))


library(lemon)
grid_arrange_shared_legend(justherb,together, ncol=2)


###
#check predictors of structure through space and time
###

both_pop_less <- both_pop %>% filter(state != "Michigan") %>% filter(state != "Indiana" ) %>% filter(env != "") #removing states with barely any obs
both_pop_less$time_span<-as.factor(as.numeric(both_pop_less$time_span))
options(scipen=0)

full<-lm(data=both_pop_less, K1 ~ long + lat + env + time_span + time_span:long + time_span:state + time_span:env:long )
Anova(full,type=3, singular.ok = F) #model results

