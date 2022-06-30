brary(dplyr)
library(ggmap)
library(scatterpie)
library(wesanderson)
library(PNWColors)
library(MetBrewer)

cg<-read.table("~/3waymerged_sampleinfo.txt",sep="\t",na.strings = c("","NA"),header=T)
head(cg)

#assign a dataset variable
cg$dataset<-"PNAS"
cg$dataset[cg$year == 2019]<-"Paired"
cg$dataset[cg$year < 2015]<-"Herbarium"


#get pop level coords
coord_byenv<-cg %>% group_by(env, long, lat, dataset) %>% dplyr:::summarise(n=n(), state=first(state))

#download map data
water<-get_stamenmap(bbox= c(left = -97, bottom = 36, 
                             right = -75, top = 45), maptype='toner-lite', zoom=6,color = "bw")

#different zoom
#water2<-get_stamenmap(bbox= c(left = -83.5, bottom = 41.8, 
#                             right = -81, top = 43), maptype='toner-background”', zoom=12)


#import color palettes
bay<-pnw_palette("Bay",7,type="continuous")
gaug<-met.brewer(name="Gauguin")

#save ggmap
p<- ggmap(water) +
  theme_nothing()

#package for multiple color schemes
library(ggnewscale)

#rename environmental variables for legend
coord_byenv$env[coord_byenv$env == "Ag"] <- "Agricultural"
coord_byenv$env[coord_byenv$env == "Nat"] <- "Natural"
cg<-cg[complete.cases(cg$env),]

#plot
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
  
#another related plot 
#Fig 1B
herb<-cg[cg$dataset == "Herbarium",]
herb_counts<-herb %>% group_by(year, env) %>% summarise(n=n())

ggplot(herb_counts, aes(year, n)) + 
  geom_col() +
  geom_rug(aes(year,y=NULL)) +
  labs(x="Collection Year", y="Sample Size") +
  theme_bw()


