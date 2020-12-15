###########Libraries#################
library(arules)
library(dplyr)
library(reshape2)
require(infotheo)
library(arulesSequences)
library(TraMineR)
library(arulesViz)
path1="D:/job/Siemens/data/"
##########FUNCTIONS###########################
Code_B_A <- function(x){
  a = as.POSIXct(x[1], format = '%Y-%m-%d %H:%M:%S')
  b = as.POSIXct(x[2], format = '%Y-%m-%d %H:%M:%S')
  if (a > b){
    return('Before')
  }else if ((a <= b) & ((a + as.numeric(x[3]) * 60) >= b)){
    return('In')
  }else if ((a + as.numeric(x[3]) * 60) < b){
    return('After')
  }
}
###################
#Load dataset
dt1 <- read.csv(paste0(path1,"/All Sites Together encoded.csv"),header =T)
dt5 <- read.csv(paste0(path1,"/Codes and Event Warning Stop classification.csv"),header =T)
#Data Prrepration
dt1['VisitStartTime'] = as.POSIXct(dt1$VisitStartTime,format = '%m/%d/%Y %H:%M')
dt1['TimeOn'] = as.POSIXct(dt1$TimeOn,format = '%m/%d/%Y %H:%M')
dt1['TimeOff'] = as.POSIXct(dt1$TimeOff,format = '%m/%d/%Y %H:%M')
dt1$difftim <- difftime(dt1$TimeOff,dt1$TimeOn,units = c("mins"))
dt1$StationID <- as.factor(dt1$StationID)
dt1$VisitId <- as.factor(dt1$VisitId)
dt1$Code <- as.factor(dt1$Code)
dt1 <- merge(dt1,dt5,by="Code")
dt1["Code_B_A"] = apply(dt1[c('VisitStartTime', 'TimeOn', 'VisitDurMinutes')], 1, Code_B_A)
dt1['visnew']=paste0(dt1$VisitId,"_",dt1$Code_B_A)
##############################
#extract association rueles baed on sequence of occurnace
trnsbf <- as(split(dt1[which(dt1$Code_B_A=="Before"),"Code"], dt1[which(dt1$Code_B_A=="Before"),"visnew"]), "transactions")
trnsaf <- as(split(dt1[which(dt1$Code_B_A=="After"),"Code"], dt1[which(dt1$Code_B_A=="After"),"visnew"]), "transactions")
trnsin <- as(split(dt1[which(dt1$Code_B_A=="In"),"Code"], dt1[which(dt1$Code_B_A=="In"),"visnew"]), "transactions")
#Apply Apriori for Before dataset
rlbf <- apriori(trnsbf, parameter = list(supp = 0.15, conf = 0.99, minlen=2))
summary(rlbf)
#Transform rules to dataset
before_1 <- as.data.frame(inspect(sort(rlbf, by = "lift")))
plot(rlbf, method="grouped",  control=list(k=50))
subrules2 <- head(sort(rlbf, by="lift"), 20)
plot(subrules2, method="graph")
#Apply Apriori for After dataset
rlaf <- apriori(trnsaf, parameter = list(supp = 0.03, conf = 0.09, minlen=3))
summary(rlaf)
#Transform rules to dataset
after_1 <- inspect(sort(rlaf, by = "lift"))
plot(rlaf, method="grouped",  control=list(k=50))
subrules2 <- head(sort(rlaf, by="lift"), 20)
plot(subrules2, method="graph")
#Apply Apriori for In dataset
rlin <- apriori(trnsin, parameter = list(supp = 0.1, conf = 0.9, minlen=2))
summary(rlin)
#Transform rules to dataset
in_rules <- inspect(head(sort(rlin, by = "lift"),400))
plot(rlin, method="grouped",  control=list(k=50))
subrules2 <- head(sort(rlin, by="lift"), 20)
plot(subrules2, method="graph")
