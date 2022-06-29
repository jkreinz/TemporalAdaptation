##########
#Figure 1
##########
cg<-read.table("~/3waymerged_sampleinfo.txt",sep="\t",na.strings = c("","NA"),header=T)
#pnas_contemp<-cg[-c(188:295),]
head(cg)
cg$dataset<-"PNAS"
cg$dataset[cg$year == 2019]<-"Paired"
cg$dataset[cg$year < 2015]<-"Herbarium"

#cg[cg$state == "Walpole",]$dataset <- "Paired"
#cg[cg$state == "Nat",]$dataset <- "Paired"

library(dplyr)
coord_byenv<-cg %>% group_by(env, long, lat, dataset) %>% dplyr:::summarise(n=n(), state=first(state))
#coord_byenv<-coord_byenv[-c(48,49,53:57),]
#coord_byenv<-coord_byenv[-c(123,124),]
#coord_byenv<-coord_byenv[-c(132),]

library(ggmap)

water<-get_stamenmap(bbox= c(left = -97, bottom = 36, 
                             right = -75, top = 45), maptype='toner-lite', zoom=6,color = "bw")

#water2<-get_stamenmap(bbox= c(left = -83.5, bottom = 41.8, 
#                             right = -81, top = 43), maptype='toner-background”', zoom=12)


library(ggmap)
library(scatterpie)
library(wesanderson)
library(PNWColors)
bay<-pnw_palette("Bay",7,type="continuous")
#devtools::install_github("BlakeRMills/MetBrewer") 
library("MetBrewer")
gaug<-met.brewer(name="Gauguin")

p<- ggmap(water) +
  #geom_point(data= owgps, aes(x = Longitude, y = Latitude, size = 2), alpha=.5, colour="darkred") +
  #geom_path(data=mdat,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #geom_path(data=ca.provinces,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #scale_color_manual(values = c(prop = "#FF0000", suseptible = "#F2AD00", nat="#00A08A")) 
  theme_nothing()

#p2<- ggmap(water2) +
  #geom_point(data= owgps, aes(x = Longitude, y = Latitude, size = 2), alpha=.5, colour="darkred") +
  #geom_path(data=mdat,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #geom_path(data=ca.provinces,aes(x=long,y=lat,group=group), colour="black",alpha=.8) +
  #scale_color_manual(values = c(prop = "#FF0000", suseptible = "#F2AD00", nat="#00A08A")) 
  #theme_nothing()

library(ggnewscale)
ggmap(water) + 
  #geom_jitter(data=coord_byenv[coord_byenv$dataset == "PNAS",],aes(y = lat, x = long, color=state), alpha=.7, size=6, height=.15, width=.15) +
  #scale_color_manual(values=c(bay[c(5:10)])) +
  #labs(color="Kreiner2019") +
  #new_scale_color() +
  geom_jitter(data=coord_byenv[coord_byenv$dataset == "Paired",],aes(y = lat, x = long, color=env), alpha=.75, size=3, height=.2, width=.15) +
  scale_color_manual(values=c(gaug[1],gaug[5])) +
  labs(colour="Contemporary \nPopulation Pairs") +
  new_scale_color() +
  geom_point(data=cg[cg$dataset == "Herbarium",],aes(y = lat, x = long, color=year, shape=env), alpha=.7, size=2.5) +
  scale_color_continuous(low = "lightgrey", high="black") +
  xlab("Longitude") +
  ylab("Latitude")
  #labs(colour="Herbarium \nSamples")


coord_byenv$env[coord_byenv$env == "Ag"] <- "Agricultural"
coord_byenv$env[coord_byenv$env == "Nat"] <- "Natural"
cg<-cg[complete.cases(cg$env),]

ggmap(water) + 

  geom_jitter(data=coord_byenv[coord_byenv$dataset == "Paired",],aes(y = lat, x = long, color=env), size=4.5, alpha=.9, height=.2, width=.15) +
  scale_color_manual(values=c("grey","grey40")) +
  labs(colour="Contemporary \nPopulation Pairs") +
  new_scale_color() +
  geom_point(data=cg[cg$dataset == "Herbarium",],aes(y = lat, x = long, color=year, shape=env), size=2.5) +
  scale_color_viridis_c() +
  xlab("Longitude") +
  ylab("Latitude")+
  labs(colour="Herbarium \nSamples")
  
herb<-cg[cg$dataset == "Herbarium",]
herb_counts<-herb %>% group_by(year, env) %>% summarise(n=n())
#herb_counts<-cg %>% group_by(year, env) %>% summarise(n=n())

ggplot(herb_counts, aes(year, n)) + 
  geom_col() +
  geom_rug(aes(year,y=NULL)) +
  labs(x="Collection Year", y="Sample Size") +
  theme_bw()


