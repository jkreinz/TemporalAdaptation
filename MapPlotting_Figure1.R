##########
#Figure 1
##########
library(dplyr)
library(ggmap)
library(scatterpie)
library(wesanderson)
library(PNWColors)
library(MetBrewer)
library(maps)
library(mapdata)
library(maptools)  #for shapefiles
library(scales)  #for transparency

cg<-read.table("~/3waymerged_sampleinfo.txt",sep="\t",na.strings = c("","NA"),header=T)
head(cg)

#assign a dataset variable
cg$dataset<-"PNAS"
cg$dataset[cg$year == 2019]<-"Paired"
cg$dataset[cg$year < 2015]<-"Herbarium"


#get pop level coords
library(dplyr)
coord_byenv<-cg %>% group_by(env, long, lat, dataset) %>% dplyr:::summarise(n=n(), state=first(state))

#get map base
library(raster)
states    <- c("New York","Pennsylvania","Maryland","West Virginia","Virginia","Kentucky","Ohio",
               "Michigan","Indiana","Illinois","Wisconsin","Minnesota","Iowa","Missouri",
               "Kansas","Nebraska","South Dakota","North Carolina","Tennessee","Mississippi", "Oklahoma",
               "Lake Michigan","Lake Ontario","Lake Superior")
provinces <- c("Ontario")

us <- getData("GADM",country="USA",level=1)
canada <- getData("GADM",country="CAN",level=1)

us.states <- us[us$NAME_1 %in% states,]
ca.provinces <- canada[canada$NAME_1 %in% provinces,]

#lakes
lakes <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

#ocean
ocean <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'ocean', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)

# rivers
rivers <- rnaturalearth::ne_download(scale = 10, 
                                     type = 'rivers_lake_centerlines', 
                                     category = 'physical')  %>% 
  sf::st_as_sf(lakes110, crs = 4269)

library("rnaturalearth")
install.packages("rnaturalearthhires")
library("rnaturalearthhires")
library(rnaturalearth)
library(dplyr)
library(raster)
library(sf)
library(tidyverse)
library(ggrepel)

#get US states outlines
in_sf <- ne_states(geounit = "United States of America",
                   returnclass = "sf")

uss <- st_as_sf(us.states) %>% 
  mutate(
    lon = map_dbl(geometry, ~st_centroid(.x)[[1]]),
    lat = map_dbl(geometry, ~st_centroid(.x)[[2]]))

#move Michigan label and add Ontario label
uss[uss$NAME_1=="Michigan",]$lon<--85.0554
uss[uss$NAME_1=="Michigan",]$lat<-44.00902
uss<-uss %>% add_row(NAME_1 = "Ontario", lon = -78.5554, lat=45)

#rename environmental variables for legend
coord_byenv$env[coord_byenv$env == "Ag"] <- "Agricultural"
coord_byenv$env[coord_byenv$env == "Nat"] <- "Natural"
cg$env[cg$env == "Ag"] <- "Agricultural"
cg$env[cg$env == "Nat"] <- "Natural"
cg$env[cg$env == "Dist"] <- "Disturbed"

#cg<-cg[complete.cases(cg$env),]

coord_byenv$env<-as.factor(coord_byenv$env)
cg$env<-as.factor(cg$env)
cg$env<-factor(cg$env, levels = c("Agricultural", "Natural", "Disturbed"))


#generate full map base
library(ggthemes)
plain1<- 
  ggplot()+
  geom_path(data=us.states,aes(x=long,y=lat,group=group))+
  geom_path(data=ca.provinces, aes(x=long,y=lat,group=group), fill="grey")+
  coord_map() +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = ocean,
          mapping = aes(geometry = geometry),
          color = "black")  +
  geom_sf(data = rivers,
          mapping = aes(geometry = geometry),
          color = "grey80",alpha=.75)  +
  theme_nothing() +
  scale_x_continuous(limits=c(-97,-75)) +
  scale_y_continuous(limits=c(36,45)) +
  coord_sf(xlim=c(-96,-75)) +
  geom_text(
    data = uss,
    aes(x = lon,
        y = lat,
        label = NAME_1),cex=3, alpha=.8) 


#import color palettes
bay<-pnw_palette("Bay",7,type="continuous")
gaug<-met.brewer(name="Gauguin")

#package for multiple color schemes
library(ggnewscale)

#rename environmental variables for legend
coord_byenv$env[coord_byenv$env == "Ag"] <- "Agricultural"
coord_byenv$env[coord_byenv$env == "Nat"] <- "Natural"
cg<-cg[complete.cases(cg$env),]

#plot map base with data points
plain1 + 
  geom_jitter(data=coord_byenv[coord_byenv$dataset == "Paired",],aes(y = lat, x = long, color=env, shape=env), size=5, alpha=.9, height=.2, width=.15) +
  scale_color_manual(values=c("grey","grey40")) +
  scale_shape_manual(values=c("circle","square","triangle")) +
  labs(colour="Contemporary\nPaired Population\nEnvironment") +
  new_scale_color() +
  geom_jitter(data=cg[cg$dataset == "Herbarium",],aes(y = lat, x = long, color=year, shape=env), size=3, height=.1,width=.05, alpha=.9) +
  scale_color_viridis_c() +
  xlab("Longitude") +
  ylab("Latitude")+
  labs(colour="Herbarium\nCollection\nYear",shape="Herbarium\nCollection\nEnvironment") 

#another related plot 
#Fig 1B
herb<-cg[cg$dataset == "Herbarium",]
herb_counts<-herb %>% group_by(year, env) %>% summarise(n=n())

p2<-ggplot(herb_counts, aes(year, n)) + 
  geom_col() +
  geom_rug(aes(year,y=NULL)) +
  labs(x="Collection Year", y="Sample Size") +
  theme(aspect.ratio = 0.5) +
  theme_bw()

row1<-plot_grid(p1, nrow=1)
row2<-plot_grid(p2, nrow=1, scale=.8)

library(cowplot)
plot_grid(row1, row2, labels = "AUTO", rel_heights=c(4,1),
          nrow=2) 


library(gridExtra)
grid.arrange(p1,p2, nrow=2, widths=c(3,1),nrow=2)
