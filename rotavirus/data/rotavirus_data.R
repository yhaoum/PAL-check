###############################################
# In this code we calculate the weekly number of new rotavirus cases in Germany
# between the years 2001-2008 stratified by age, namely age groups 0-4, 5-59 and 60+, 
# scaled up the the underreporting rates as inferred in Weidemann et al. 2013. 
# The data loaded is available at the following Github repository:
# https://github.com/weidemannf/thesis_code/tree/master/Mode_inference_paperI
# Author: Theresa Stocks
# Date: June 28, 2017
#################################################




# init
rm(list = ls())
#setwd("~/Dropbox/two-age-strata/Felix_data/thesis_code-master/Mode_inference_paperI")

## [[]] used to cite columns of a data.frame
EFSdata <- function(n){
  EFSdat<-numeric(0)
  for (i in 1:n){
    direct<-paste("source/RotadataEFS200",i,".txt")
    rotadata <- read.table(file=direct)
    for (j in 1:min(52,length(rotadata))){EFSdat<-cbind(EFSdat,rotadata[[j]])}
  }
  return(EFSdat)
}

# number of age groups 1:5 children, 6:8 adults, 9:10 elderly 
# Actually, year = 8
nrow(EFSdata(8))


# we consider data form 2001-2008
year <- 8
child <- seq(1,5,by=1)
adult <- seq(6,8,by=1)
elderly <- seq(9,10, by=1)

EFS <- matrix(nrow=3,ncol= ncol(EFSdata(8)))
EFS[1,] <- colSums(EFSdata(year)[child,])
EFS[2,]<- colSums(EFSdata(year)[adult,])
EFS[3,] <- colSums(EFSdata(year)[elderly,])
EFS

WFSdata <- function(n){
  WFSdat<-numeric(0)
  for (i in 1:n){
    direct<-paste("source/RotadataWFS200",i,".txt")
    rotadata <- read.table(file=direct)
    for (j in 1:min(52,length(rotadata))){WFSdat<-cbind(WFSdat,rotadata[[j]])}
  }
  return(WFSdat)
}

WFS <- matrix(nrow=3,ncol= ncol(EFSdata(8)))
WFS[1,] <- colSums(WFSdata(year)[child,])
WFS[2,] <- colSums(WFSdata(year)[adult,])
WFS[3,] <- colSums(WFSdata(year)[elderly,])
WFS

# change of reporting behaviour in beginning of 2005, so from beg 2001- end 2004 
# constant and from beg 2005- end 2008
years_till_change <- 4
time_unit <- 52
time_change <- years_till_change* time_unit
before <- seq(1,time_change,by=1)
after <- seq(time_change+1, ncol(EFSdata(8)), by=1)

# averaged posterior density of fitted underreportings rates from Weidemann et al. 2013 before and after 2005
EFS_rep_before <- 0.19 # Before 2005 underreporting rate
EFS_rep_after <- 0.241 # After 2005 underreporting rate
WFS_rep_before <- 0.043 # Before is lower than After
WFS_rep_after <- 0.063
#######################################################################

# Stocks' data is (EFS_{0:208}/0.19,EFS_{209:416}/0.241) + 
#              (WFS_{0:208}/0.043,WFS_{209:416}/0.063)
EFS_no_rep <- round(cbind(EFS[, before]/EFS_rep_before, EFS[, after]/EFS_rep_after))
WFS_no_rep <- round(cbind(WFS[, before]/WFS_rep_before, WFS[, after]/WFS_rep_after))

rotadata_mine <- t(EFS_no_rep + WFS_no_rep)
# read-in Stocks' data
Stocks <- read.table("rotavirus.txt")[, -1]
# Check
all(rotadata_mine == Stocks) #True

# PAL data is EFS + WFS
rota_pal <- t(EFS + WFS)
# read-in PAL's data
load("real_rotavirus_metadata.Rdata")
# Check
all(rota_pal == realdat) #True