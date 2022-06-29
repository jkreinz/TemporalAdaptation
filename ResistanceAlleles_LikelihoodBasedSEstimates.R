#Likelihood estimator of selection given two lineages
#Feb 2022
#Assumes an input file with {date, lineage, n}, 
#where n stands for the observed number of sequences of that lineage on that date.

#load libraries:
library(tidyverse)    #data wrangling and plotting
#install.packages("bbmle")
library(bbmle)

#set directory:
#setwd("/Users/otto/Documents/Students/JuliaKreiner/SelectionMigration/Data") #<- CHANGE! (MACs I think use \ ) to working directory

#load lineage data with {date, Pango.lineage, n}
mydata <-  read.csv("seven_herbicideresistant_loci_genotypes_throughtime.csv")
names(mydata)[6]<-"time"
mydata <- filter(mydata, mydata$variable == "ALS653")

#Set a starting date and a final date:
#Note that the startdate shouldn't be too much before both alleles become common
#or rare migration events that die off could throw off the estimation procedure 
#(so that the parameter estimates account for the presence of those alleles long in the past).
startdate<-"1960"
lastdate<-"2018"

#filter data to after that starting date
mydata <- filter(mydata, mydata$time >= 1960)

#turn genotypic counts into the number of alleles of type 1 and type 2 (for use in existing code)
data1 <- filter(mydata, mydata$value == 0)
data1$value <- 2 #homozygotes for the sensitive allele
data2 <- filter(mydata, mydata$value == 1)
data2$value <- 2 #homozygotes for the resistant allele
datahet <- filter(mydata, mydata$value == 0.5)
datahet$value <- 1
data1 <- rbind(data1,datahet)
data2 <- rbind(data2,datahet)

data1$time <- data1$time - 2018
data2$time <- data2$time - 2018


#join lists in a dataframe to plot proportions
timestart <- 1960-2018
timeend <- 0
toplot <- data.frame(time = seq.int(timestart,timeend))
toplot$n1 <- data1$value[match(toplot$time,data1$time)]
toplot$n2 <- data2$value[match(toplot$time,data2$time)]
toplot[is.na(toplot)] <- 0 #All NA's are individuals sampled with the other alleles and count as 0

#NOT USED - end date is a good place to have t=0 for the herbicide resistance work.
#To aid in the ML search, we rescale time to be centered as close as possible
#to the midpoint (p=0.5), to make sure that the alleles are all segregating at the reference date.
#If we have t=0 at when p is near 0 or 1, then the likelihood surface is very flat.
#refdate<-which(abs(toplot$n2/(toplot$n1+toplot$n2)-0.5)==min(abs(toplot$n2/(toplot$n1+toplot$n2)-0.5),na.rm=TRUE))
#timeend <- (timeend-timestart)-refdate
#timestart <- -refdate
#toplot$time <- seq.int(timestart,timeend)
#data1$time <- data1$time + (timeend-timestart)-refdate
#data2$time <- data2$time + (timeend-timestart)-refdate

plot(y=toplot$n2/(toplot$n1+toplot$n2),x=toplot$time,xlab="Time",ylab="proportion",ylim=c(0,1))

################################
# Using mle2 and profile in BBMLE
################################
#It looks like the BBMLE package performs well and gives
#confidence intervals for the parameters.  Here, we have to flip the sign
#of the log-likelihood directly for use with mle2 (can't send control=list(fnscale=-1) through?).
binfunc <- function(p,s){
  -(sum(data1$value*log((1-p)/((1-p)+p*exp(s*data1$time))))+
      sum(data2$value*log(p*exp(s*data2$time)/((1-p)+p*exp(s*data2$time)))))
}

startpar<-list(p=0.5, s=0.1)
bbml<-mle2(binfunc, start = startpar)
bbml

#Confidence intervals (uniroot should be most accurate, but it doesn't work for the herbicide data)  
#confint(bbml) # based on inverting a spline fit to the profile 
#
#confint(bbml,method="quad") # based on the quadratic approximation at the maximum likelihood estimate
#
#confint(bbml,method="uniroot") # based on root-finding to find the exact point where the profile crosses the critical level
#
#Interesting way of profiling the likelihood graphically (if the above worked!)
bbprofile<-profile(bbml)
plot(bbprofile)
proffun(bbml)

################################
# Manually determining the CI
################################
# Since the BBMLE package failed to find the confidence interval by
# profile likelihood, we can do it manually
#Make usre that the sequence of s values used is broad enough
#(expand if a boundary value is obtained)
forceprofile <- data.frame(s = seq(-0.1,3,0.001))
startpar<-list(p=0.5)
counter<-NULL
 for (j in 1:length(forceprofile$s))  {
  counter<-counter+1
  s<-forceprofile$s[j]
  binfuncfixs <- function(p){
    -(sum(data1$value*log((1-p)/((1-p)+p*exp(s*data1$time))))+
        sum(data2$value*log(p*exp(s*data2$time)/((1-p)+p*exp(s*data2$time)))))
  }
  bbml<-mle2(binfuncfixs, start = startpar)
  forceprofile$ml[j]<-bbml@min
}

minml<-min(forceprofile$ml)
temp<-filter(forceprofile, forceprofile$ml < minml + qchisq(0.95, 1)/2)
min(temp$s)
max(temp$s)


