library(data.table)
library(ggplot2)

anc_files <- list.files(pattern="*ancestry_scaf*", full.names=TRUE, recursive=T)
fam_files <- list.files(path="/ohta1/julia.kreiner/waterhemp/commongarden/regional_vcfs/bychrom/plink/lamp_analysis/contemp_fam",pattern="*fam", full.names=TRUE, recursive=T)
fam_files_rep<-rep(fam_files,16)


datalist = list()
#perform stats by sample and scaffold from LAMP output
for (i in 1:length(fam_files_rep)) {
	
	fam<-fread(fam_files_rep[i])		
	anc<-fread(anc_files[i])
	toDelete <- seq(0, nrow(anc), 2)
	anc2 <-  anc[-toDelete, -1]
	
	anc_prop<-rowSums(anc2,na.rm=T)/ncol(anc2)
	pop<-fam$V1
	samp<-fam$V2
	scaf<-strsplit(anc_files, "_")[[i]][3]

	data<-data.frame(pop=pop, samp=samp, anc_prop=anc_prop, scaf=scaf)	
	datalist[[i]]<-data
}
#merge across samples/scafs
big_data = do.call(rbind, datalist)

#get stats by pop, sample, and scaffold
bysamp <- big_data %>% group_by(samp) %>% dplyr::summarise(mean_anc=mean(anc_prop), pop=first(pop), lowerCI=quantile(anc_prop,probs=0.05), upperCI=quantile(anc_prop,probs=0.95))

bypop <- big_data %>% group_by(pop,scaf) %>% dplyr::summarise(mean_anc=mean(anc_prop), pop=first(pop), lowerCI=quantile(anc_prop,probs=0.05), upperCI=quantile(anc_prop,probs=0.95))

byscaf <-big_data %>% group_by(scaf) %>% dplyr::summarise(mean_anc=mean(anc_prop), pop=first(pop), lowerCI=quantile(anc_prop,probs=0.05), upperCI=quantile(anc_prop,probs=0.95))


write.table(bysamp, "mean_anc_bysample_avgacrossindsandscaf.txt",quote=F,row.names=F,col.names=F)
write.table(bypop, "mean_anc_bypop_byscaf_avgacrossindsandscaf.txt",quote=F,row.names=F,col.names=F)
write.table(byscaf, "mean_anc_byscaf_avgacrossindsandscaf.txt",quote=F,row.names=F,col.names=F)


#read in positional info to window
pos_files <- list.files(pattern="*justpos", full.names=TRUE, recursive=T)
pos_files_rep <- rep(pos_files, each=34)

datalist2 = list()
#perform stats by sample and scaffold from LAMP output
for (i in 1:length(fam_files_rep)) {

        fam<-fread(fam_files_rep[i])
        anc<-fread(anc_files[i])
	pos<-fread(pos_files_rep[i])

        toDelete <- seq(0, nrow(anc), 2)
        anc2 <-  anc[-toDelete, -1]
	anc2 <- anc2[,1:((ncol(anc2)-1))]
	print(row.names(anc2))
	anc3 <- as.data.frame(t(anc2))
	names(anc3)<-fam$V2
	anc3$pos<-pos$V1	
	anc_windowed <- anc3 %>% mutate(win=floor(pos/100000)) %>% group_by(win) %>% dplyr::summarise_at(which(sapply(anc3, is.numeric)), mean)
	anc_windowed$scaf<-strsplit(anc_files, "_")[[i]][3]
	pos_win<- pos %>%  mutate(win=floor(V1/100000)) %>% group_by(win) %>% dplyr:::summarise(n_win=n()) %>% filter(n_win > 20)

	datalist2[[i]]<-semi_join(anc_windowed, pos_win, by="win")
}

library(dplyr) 
library(purrr)

all_scafs<-list()
for (i in 1:16) {
print(paste("scaf",i,sep=""))
cond <- sapply(datalist2, function(df) all(df$scaf == paste("scaf",i,sep="")))
scaf1 <- datalist2[cond]
big_data2 = do.call(cbind, scaf1)
test<-big_data2[!duplicated(as.list(big_data2))]
all_scafs[[i]]<-melt(test, id.vars = c("win","scaf"))

}


ancwinds_allscafs = do.call(rbind, all_scafs)
write.table(ancwinds_allscafs, "ancwindows_byscaf_bysample_20snpthresh.txt", col.names=T, row.names=F,quote=F)



pos_files <- list.files(pattern="*justpos", full.names=TRUE, recursive=T)

for (i in 1:length(pos_files)) {

        pos<-fread(pos_files[i])
	pos_win<- pos %>%  mutate(win=floor(V1/100000)) %>% group_by(win) %>% dplyr:::summarise(n_win=n())
	













