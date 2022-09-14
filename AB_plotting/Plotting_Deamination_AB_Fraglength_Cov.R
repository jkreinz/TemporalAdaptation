#compiled output from MapDamage of deamination levels
ct<-read.table("CtoT.all")

library(ggplot2)
meta<-read.csv("Downloads/Curated Herbarium Tissue Samples - Sheet12.csv")
names(ct)<-c("readpos","CtoT","Sample")
ct$Sample<-gsub("HB0","HBO",ct$Sample)
ct_meta<-inner_join(meta,ct,by="Sample")

v1<-ggplot(data=ct_meta, aes(readpos, CtoT, color=Collection.Year, group=Sample)) +
  geom_point(alpha=.5, cex=2) +
  geom_line(alpha=.3) +
  ylab("Percent C to T\ndeamination") +
  xlab("Read Position from 5p") +
  theme_bw() +
  scale_color_viridis_c()

ga<-read.table("GtoA.all")
ggplot(data=ga, aes(V1, V2, color=V3)) +
  geom_point(alpha=.5, cex=2) +
  geom_line(alpha=.3) +
  ylab("Percent G to A Deamination") +
  xlab("Read Position from 3p") +
  theme_bw() +
  theme(legend.position = "none")   

library(tidyverse)


v2<-ct_meta %>% filter(readpos==1) %>%
  ggplot(aes(Collection.Year, CtoT, color=Collection.Year)) +                              
  ylab("Percent C to T\ndeamination at first base") +
  geom_smooth(method="lm") +
  geom_point(alpha=.8, cex=2) + 
  xlab("Year") +
  theme_bw() +
  scale_color_viridis_c()


v3<-ct_meta %>% filter(readpos==1) %>%
  ggplot(aes(Collection.Year, FinalCov, color=Collection.Year)) +                              
  ylab("Coverage") +
  geom_smooth(method="lm") +
  geom_point(alpha=.8, cex=2) + 
  xlab("Year") +
  theme_bw()  +  scale_color_viridis_c()

frag<-read.csv("Downloads/Curated Herbarium Tissue Samples - Sheet15.csv")
names(frag)[1]<-"Sample"
frag_meta<-inner_join(meta,frag,by="Sample")

v4<-frag_meta %>% 
  ggplot(aes(Collection.Year, fragmentlength, color=Collection.Year)) +                              
  ylab("Fragment Length") +
  geom_smooth(method="lm") +
  geom_point(alpha=.8, cex=2) + 
  xlab("Year") +
  theme_bw()  +  scale_color_viridis_c()

lm1<-lm(data=frag_meta, FinalCov~Collection.Year)
summary(lm1)
anova(lm1)

library(gridExtra)
library(lemon)
grid_arrange_shared_legend(v1,v2,v3,v4, ncol=2 , nrow=2,position = "bottom")

AB_head<-read.table("AB_headerorder.txt")
AB_alt<-read.table("alternate_AB_bygeno.txt")
AB_ref<-read.table("reference_AB_bygeno.txt")
AB_anc<-read.table("herbarium_ancestyrinformative_ABbygeno.txt")ymu

AB_all<-rbind(AB_alt,AB_ref)
AB_focalmean<-colMeans(AB_all[,-c(1:4)],na.rm = T)
AB_focalse<-sapply(AB_all[,-c(1:4)],function(x)sd(x,na.rm = T)/sqrt(length(na.omit(x))))

AB_ancmean<-colMeans(AB_anc[,-c(1:4)],na.rm = T)
AB_ancse<-sapply(AB_anc[,-c(1:4)],function(x)sd(x,na.rm = T)/sqrt(length(na.omit(x))))

df2<-data.frame("Sample"=AB_head, "AB_focal"=AB_focalmean, "AB_anc"=AB_ancmean,
                "AB_focal_se"=AB_focalse,
                "AB_anc_se"=AB_ancse)
names(df2)[1]<-"Sample"
df2$Sample<-gsub("HB0","HBO",df2$Sample)
AB_meta<-inner_join(df2,meta,by="Sample")

AB_meta2<-AB_meta[complete.cases(AB_meta[,c(1:5)]),]

v5<-ggplot(data=AB_meta, aes(Collection.Year, AB_focal, color=Collection.Year)) +
  geom_point() +
  geom_smooth(method="lm") +
  scale_color_viridis_c() +
  theme_bw() +
  ylim(0.3,.7) +
  geom_errorbar(ymin=AB_meta$AB_focal-AB_meta$AB_focal_se,ymax=AB_meta$AB_focal+AB_meta$AB_focal_se) +
  labs(y="Mean Allelic Bias by Sample\n(Focal Ag-Alleles)")

v6<-ggplot(data=AB_meta, aes(Collection.Year, AB_anc, color=Collection.Year)) +
  geom_point() +
  geom_smooth(method="lm") +
  scale_color_viridis_c() +
  theme_bw() +
  ylim(0.3,0.7) +
  geom_errorbar(ymin=AB_meta$AB_anc-AB_meta$AB_anc_se,ymax=AB_meta$AB_anc+AB_meta$AB_anc_se) +
  labs(y="Mean Allelic Bias by Sample\n(Ancestry Informative Alleles)")

summary(lm1<-lm(data=AB_meta, AB_focal ~ Collection.Year))
summary(lm1<-lm(data=AB_meta, AB_anc ~ Collection.Year))
cor(AB_meta$AB_anc, AB_meta$Collection.Year)
cor(AB_meta$AB_focal, AB_meta$Collection.Year, use="complete.obs")

library(cowplot)
# extract the legend from one of the plots
legend_b <- get_legend(
  v2 + 
    guides(fill = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

prow<-plot_grid(
  v1 + theme(legend.position="none"),
  v2 + theme(legend.position="none"),
  v3 + theme(legend.position="none"),
  v4 + theme(legend.position="none"),
  v5 + theme(legend.position="none"),
  v6 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C","D","E","F"),
  #hjust = -1,
  nrow = 3,
  ncol = 2,
  scale=.9
)


plot_grid(prow, legend_b, nrow=2,rel_heights = c(1, .05))
