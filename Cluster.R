###########Libraries#################
require(clusterCrit)
library(ggplot2)
library(reshape2)
library(car)
require(infotheo)
library(statnet)
library(GGally)
library(igraph)
path1="D:/job/Siemens/data/"
##############################################
##########FUNCTIONS###########################
##PLot of clustering (save the plots in WD)
#cln is nu,ber of clusters which has been found, name1 is file name, rec.vis contains visits and associated clusters
plt1 <- function(dt1,rec.vis,name1,cln){
  dtvis <- merge(dt1,rec.vis,by="VisitId")
  tb1 <- as.data.frame(table(dtvis$Park_Name,dtvis$VisCluster))
  tb2 <- as.data.frame.matrix(table(dtvis$Park_Name,dtvis$VisCluster))
  rs1 <- rowSums(as.matrix(tb2))
  tb1$row1 <- tb1$Freq/rep(rs1,cln)
  rs2 <- colSums(as.matrix(tb2))
  rep3=c()
  for (i in cln:1){
    obj1=rep(rs2[i],37)
    rep3=  append(obj1,rep3)
  }
  tb1$col1 <- tb1$Freq/rep3
  tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
  p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=row1,size=row1),colour = "Red")+scale_size_continuous(range = c(0.5,8))+
    labs(x="Clusters",y="Parks",size="Frequency",alpha="Frequency") +
    ggtitle("Bubble plot of visit patterns(clusters) - Normalized by Park")+
    theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
          plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
          axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
          axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
          legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
          legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
  ggsave(paste0(name1,"_NorP.jpg"),width = 10, height = 10)
  p1 <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Red")+scale_size_continuous(range = c(0.5,8))+
    labs(x="Clusters",y="Parks",size="Frequency",alpha="Frequency") +
    ggtitle("Bubble plot of visit patterns(clusters) - Normalized by cluster")+
    theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
          plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
          axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
          axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
          legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
          legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
  ggsave(paste0(name1,"_NorC.jpg"),width = 10, height = 10)
}
##Clustering Function
kmcl <- function(data, kmin, kmax, rnd = 1357)
{
  DB <- c()
  for (i in 1:ncol(data))
  {
    data[, i] <- as.numeric(data[, i])
  }
  for (j in kmin:kmax)
  {
    set.seed(rnd)
    km <- kmeans(data, centers = j)
    # km <- kmeans(data,centers=i,nstart = 10)
    data1 <- as.matrix(data)
    # Computing the Davies-Bouldin
    DB[j] <- intCriteria(data1, km$cluster, crit = c("Davies_Bouldin"))
    print("K is equal to=")
    print(j)
  }
  return(DB)
}
## Davies-Bouldin Plot
dbplot <- function(obj, tit1 = "Null")
{
  plot(2:15, unlist(obj), type = "o", col = "black", ylim = c(0, 3), 
       lwd = 2, ylab = "DB Vlaue", xlab = "Number of Clusters", cex.lab = 1.5, 
       cex.axis = 0.75)
  grid(nx = 25, ny = 25, col = "lightgray", lty = "dotted", lwd = par("lwd"), 
       equilogs = TRUE)
  axis(side = 1, at = seq(2, 15, by = 1), cex.axis = 0.75)
  box(col = "black", lwd = 2)
  title(tit1)
}
# PCA Function
PCfun <- function(dt3, pcnum = 9)
{
  for (i in 1:ncol(dt3))
  {
    dt3[, i] <- as.numeric(dt3[, i])
  }
  set.seed(1357)
  pc1 <- prcomp(dt3)
  plt <- plot(pc1$sdev, type = "l", ylab = "Standard Deviation of PCAs", 
              xlab = "Index", lwd = 2)
  box(col = "black", lwd = 2)
  return(list(dt = as.matrix(pc1$x[, 1:pcnum]), plt, devs = pc1$sdev))
}
# Create Code_B_A varaible to find out weather the code happened
# after,before or within visit
Code_B_A <- function(x)
{
  a <- as.POSIXct(x[1], format = "%Y-%m-%d %H:%M:%S")
  b <- as.POSIXct(x[2], format = "%Y-%m-%d %H:%M:%S")
  if (a > b)
  {
    return("Before")
  } else if ((a <= b) & ((a + as.numeric(x[3]) * 60) >= b))
  {
    return("In")
  } else if ((a + as.numeric(x[3]) * 60) < b)
  {
    return("After")
  }
}
# Normalization
norm1 <- function(x)
{
  (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}
##########Import dataset#######################
dt1 <- read.csv(paste0(path1,"/All Sites Together encoded.csv"),header =T)
dt5 <- read.csv(paste0(path1,"/Codes and Event Warning Stop classification.csv"),header =T)
dt6 <- read.csv(paste0(path1,"/Kcode.csv"),header =T)
#Data Prrepration
dt1$VisitStartTime <- as.POSIXct(dt1$VisitStartTime, format = "%m/%d/%Y %H:%M")
dt1$TimeOn <- as.POSIXct(dt1$TimeOn, format = "%m/%d/%Y %H:%M")
dt1$TimeOff <- as.POSIXct(dt1$TimeOff, format = "%m/%d/%Y %H:%M")
dt1$difftim <- difftime(dt1$TimeOff, dt1$TimeOn, units = c("mins"))
dt1$StationID <- as.factor(dt1$StationID)
dt1$VisitId <- as.factor(dt1$VisitId)
dt1$Code <- as.factor(dt1$Code)
dt1 <- merge(dt1, dt5, by = "Code")
dt1["Code_B_A"] <- apply(dt1[c("VisitStartTime", "TimeOn", "VisitDurMinutes")],1, Code_B_A)
dt1 <- merge(dt1, dt6, by = "Code")
###PART 1
##############################################
###clustering VisitIds########################
###In clustering section we could creat functions but number of parameters were large
###We did not provide trial and error part of our experiment here
##Aggregate based on the StationId and VisitId for all the error codes
##Records are VisitId and column are codes frequency
#Tabulate and discretize dataset
mlt1 = melt(dt1,id.vars =  names(dt1)[c(2:21)])
dt2 <- dcast(mlt1, VisitId  ~ value, value.var = c("VisitId"),fun.aggregate = length)
dt3.std <- dt2[, -1]
#Standardizing the data
for (i in 1:ncol(dt3.std))
{
  dt3.std[, i] <- infotheo::discretize(X = dt3.std[, i], disc = "equalfreq", 
                                       nbins = 5)
}
#Apply PCA on dataset
pcadt.sd <- PCfun(dt3.std, pcnum = 8)
#Apply K-means
set.seed(1357)
km1 <- kmeans(pcadt.sd$dt, centers = 9)
table(km1$cluster)
#Craete the dataset contains visit and associated clusters ans save it
rec.vis <- cbind(dt2, km1$cluster)[, c(1, 644)]
names(rec.vis) <- c("VisitId", "VisCluster")
#######Unitil here is required for PART 11
write.csv(rec.vis, "recvis.csv", row.names = FALSE)
#Create a unique data from parks, visitid, statitionID and visitType
vars1 <- c("Park_Name", "VisitId", "StationID", "VisitType")
fn1 <- unique(dt1[, which(names(dt1) %in% vars1)])
fn2 <- merge(fn1, rec.vis, by = "VisitId")
#Export graphs of parks and clsuters in which frequency is visitid
plt1(fn1, rec.vis, name1 = "All_VisitIds_Visitfreq", cln = 9)
#Export graphs of parks and clsuters in which frequency is code frequency in visitid
plt1(dt1, rec.vis, name1 = "All_VisitIds_codefreq", cln = 9)
#Some cross tables
table(fn2$VisCluster, fn2$VisitType)
#Graph for visit Type vs clusters and frequency is visitIds
fn1 <- unique(dt1[, which(names(dt1) %in% vars1)])
fn1 <- merge(fn1, rec.vis, by = "VisitId")
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster, fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster, fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1, 4)
rs2 <- colSums(as.matrix(tb2))
rep3 <- c()
for (i in 4:1)
{
  obj1 <- rep(rs2[i], 9)
  rep3 <- append(obj1, rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc = "equalfreq", nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by Visit Type")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_All data _visitfreq.jpg",width = 8, height = 8)
#Graph for visit Type vs clusters and frequency is frequency of codes codes in visitIds
fn1 <- merge(dt1,rec.vis,by="VisitId")
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster,fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster,fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,4)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 4:1){
  obj1=rep(rs2[i],9)
  rep3= append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by Visit Type")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_All data_codefreq.jpg",width = 8, height = 8)
###PART 2
##############################################
###clustering VisitIds########################
##Records are VisitId (filtered by codes happened before visit) and column are codes frequency
#Tabulate and discretize dataset
code_bai <- dt1[which(dt1$Code_B_A=="Before"),]
mlt1 = melt(code_bai,id.vars =  names(dt1)[c(2:21)])
dt2 <- dcast(mlt1, VisitId  ~ value, value.var = c("VisitId"),fun.aggregate = length)
#Standardizing the data
dt3.std <- dt2[,-1]
for (i in 1:ncol(dt3.std)){
  dt3.std[,i] <- infotheo::discretize(X = dt3.std[,i], disc="equalfreq",nbins = 5)
}
#Apply PCA on dataset
pcadt.sd <- PCfun(dt3.std,pcnum = 9)
#Apply K-means
set.seed(1357)
km1 <- kmeans(pcadt.sd$dt,centers = 7)
table(km1$cluster)
rec.vis <- cbind(dt2,km1$cluster)[,c(1,ncol(dt2)+1)]
names(rec.vis) <- c("VisitId","VisCluster")
#Craete the dataset contains visits and associated clusters and save it
write.csv(rec.vis,"recvis_Before.csv",row.names = FALSE)
#Create a unique data from parks, visitid, statitionID and visitType
vars1 <- c("Park_Name","VisitId","StationID","VisitType")
fn1 <- unique(dt1[which(dt1$Code_B_A=="Before"),which(names(dt1) %in% vars1)])
#Export graphs of parks and clsuters in which frequency is visitid
plt1(fn1 ,rec.vis ,name1 = "Before_VisitIds_visitfreq",cln = 7)
#Export graphs of parks and clsuters in which frequency is code frequency in visitid
plt1(code_bai ,rec.vis ,name1 = "Before_VisitIds_codefreq",cln = 7)
#Some cross tables
fn1 <- merge(fn1, rec.vis, by = "VisitId")
table(fn1$VisCluster,fn1$VisitType)
#Graph for visit Type vs clusters and frequency is visitIds
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster,fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster,fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,4)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 4:1){
  obj1=rep(rs2[i],7)
  rep3= append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by cluster")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_visitfreq_Before.jpg",width = 8, height = 8)
#Graph for visit Type vs clusters and frequency is frequency of codes codes in visitIds
fn1 <- merge(code_bai,rec.vis,by="VisitId")
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster,fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster,fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,4)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 4:1){
  obj1=rep(rs2[i],7)
  rep3= append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by cluster")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_codefreq_Before.jpg",width = 8, height = 8)
#Plot of visit Id clusters vs stations in a specific park (based on frequency visitid)
#You have to change the park number in prk1 to get different park's plot
prk1=19
nm1=as.character(sort(unique(dt1$Park_Name))[prk1])
fn1 <- unique(dt1[which(dt1$Code_B_A=="Before"),which(names(dt1) %in% vars1)])
fn1 <- merge(fn1, rec.vis, by = "VisitId")
kn2 = fn1[fn1$Park_Name == nm1,]
kn2$StationID = as.numeric(kn2$StationID)
kn3 <- as.data.frame(table(kn2$StationID,kn2$VisCluster))
names(kn3)[3] <- "Frequency"
p <- ggplot(kn3, aes(Var1,Var2))
p+geom_point(aes(alpha=Frequency,size=Frequency),colour = "Red")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Stations",y="Pattern (Cluster)") +
  ggtitle(paste0("Bubble plot of visit clusters frequency in ",nm1))+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 8,angle = 90,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave(paste0(nm1,"_plot of clusters and statitions_Before.jpg"),width = 8, height = 8)
###PART 3
##############################################
###clustering VisitIds########################
##Records are VisitId (filtered by codes happened after visit) and column are codes frequency
#Tabulate and discretize dataset
code_bai <- dt1[which(dt1$Code_B_A=="After"),]
mlt1 = melt(code_bai,id.vars =  names(dt1)[c(2:21)])
dt2 <- dcast(mlt1, VisitId  ~ value, value.var = c("VisitId"),fun.aggregate = length)
dt3.std <- dt2[,-1]
#Standardizing the data
for (i in 1:ncol(dt3.std)){
  dt3.std[,i] <- infotheo::discretize(X = dt3.std[,i], disc="equalfreq",nbins = 5)
}
#Apply PCA on dataset
pcadt.sd <- PCfun(dt3.std,pcnum = 9)
#Apply K-means
set.seed(1357)
km1 <- kmeans(pcadt.sd$dt,centers = 9)
table(km1$cluster)
rec.vis <- cbind(dt2,km1$cluster)[,c(1,ncol(dt2)+1)]
#Craete the dataset contains visits and associated clusters ans save it
names(rec.vis) <- c("VisitId","VisCluster")
write.csv(rec.vis,"recvis_After.csv",row.names = FALSE)
#Export graphs of parks and clsuters in which frequency is visitid
vars1 <- c("Park_Name","VisitId","StationID","VisitType")
fn1 <- unique(dt1[which(dt1$Code_B_A=="After"),which(names(dt1) %in% vars1)])
plt1(fn1 ,rec.vis ,name1 = "After_VisitIds_visitfreq",cln = 9)
#Export graphs of parks and clsuters in which frequency is code frequency in visitid
plt1(code_bai ,rec.vis ,name1 = "After_VisitIds_codefreq",cln = 9)
#Some cross tables
fn1 <- merge(fn1, rec.vis, by = "VisitId")
table(fn1$VisCluster,fn1$VisitType)
#Graph for visit Type vs clusters and frequency is visitIds
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster,fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster,fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,4)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 4:1){
  obj1=rep(rs2[i],9)
  rep3= append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by cluster")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_visitfreq_After.jpg",width = 8, height = 8)
#Graph for visit Type vs clusters and frequency is frequency of codes codes in visitIds
fn1 <- merge(code_bai,rec.vis,by="VisitId")
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster,fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster,fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,4)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 4:1){
  obj1=rep(rs2[i],9)
  rep3= append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by cluster")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_codefreq_After.jpg",width = 8, height = 8)
#Plot of visit Id clusters vs stations in a specific park (based on frequency visitid)
#You have to change the park number in prk1 to get different park's plot
prk1=19
nm1=as.character(sort(unique(dt1$Park_Name))[prk1])
fn1 <- unique(dt1[which(dt1$Code_B_A=="After"),which(names(dt1) %in% vars1)])
fn1 <- merge(fn1, rec.vis, by = "VisitId")
kn2 = fn1[fn1$Park_Name == nm1,]
kn2$StationID = as.numeric(kn2$StationID)
kn3 <- as.data.frame(table(kn2$StationID,kn2$VisCluster))
names(kn3)[3] <- "Frequency"
p <- ggplot(kn3, aes(Var1,Var2))
p+geom_point(aes(alpha=Frequency,size=Frequency),colour = "Red")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Stations",y="Pattern (Cluster)") +
  ggtitle(paste0("Bubble plot of visit clusters frequency in ",nm1))+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 8,angle = 90,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave(paste0(nm1,"_plot of clusters and statitions_After.jpg"),width = 8, height = 8)
###PART 4
##############################################
###clustering VisitIds########################
##Records are VisitId (filtered by codes happened within visit) and column are codes frequency
#Tabulate and discretize dataset
code_bai <- dt1[which(dt1$Code_B_A=="In"),]
mlt1 = melt(code_bai,id.vars =  names(dt1)[c(2:21)])
dt2 <- dcast(mlt1, VisitId  ~ value, value.var = c("VisitId"),fun.aggregate = length)
dt3.std <- dt2[,-1]
#Standardizing the data
for (i in 1:ncol(dt3.std)){
  dt3.std[,i] <- infotheo::discretize(X = dt3.std[,i], disc="equalfreq",nbins = 5)
}
#Apply PCA on dataset
pcadt.sd <- PCfun(dt3.std,pcnum = 9)
#Apply K-means
set.seed(1357)
km1 <- kmeans(pcadt.sd$dt,centers = 7)
#Craete the dataset contains visits and associated clusters ans save it
table(km1$cluster)
rec.vis <- cbind(dt2,km1$cluster)[,c(1,ncol(dt2)+1)]
names(rec.vis) <- c("VisitId","VisCluster")
write.csv(rec.vis,"recvis_In.csv",row.names = FALSE)
#Create a unique data from parks, visitid, statitionID and visitType
vars1 <- c("Park_Name","VisitId","StationID","VisitType")
fn1 <- unique(dt1[which(dt1$Code_B_A=="In"),which(names(dt1) %in% vars1)])
#Export graphs of parks and clsuters in which frequency is visitid
plt1(fn1 ,rec.vis ,name1 = "In_VisitIds_visitfreq",cln = 7)
#Export graphs of parks and clsuters in which frequency is code frequency in visitid
plt1(code_bai ,rec.vis ,name1 = "In_VisitIds_codefreq",cln = 7)
##Some cross tables
fn1 <- merge(fn1, rec.vis, by = "VisitId")
table(fn1$VisCluster,fn1$VisitType)
#Graph for visit Type vs clusters and frequency is visitIds
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster,fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster,fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,4)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 4:1){
  obj1=rep(rs2[i],7)
  rep3= append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by cluster")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_visitfreq_In.jpg",width = 8, height = 8)
#Graph for visit Type vs clusters and frequency is frequency of codes codes in visitIds
fn1 <- merge(code_bai,rec.vis,by="VisitId")
#Creating graph
tb1 <- as.data.frame(table(fn1$VisCluster,fn1$VisitType))
tb2 <- as.data.frame.matrix(table(fn1$VisCluster,fn1$VisitType))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,4)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 4:1){
  obj1=rep(rs2[i],7)
  rep3= append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=col1,size=col1),colour = "Blue")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Visit Types",y="clusters",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of VisitType and VisitId clusters - Normalized by cluster")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("VisitType and VisitId clusters_codefreq_In.jpg",width = 8, height = 8)
#Plot of visit Id clusters vs stations in a specific park (based on frequency visitid)
#You have to change the park number in prk1 to get different park's plot
prk1=19
nm1=as.character(sort(unique(dt1$Park_Name))[prk1])
fn1 <- unique(dt1[which(dt1$Code_B_A=="In"),which(names(dt1) %in% vars1)])
fn1 <- merge(fn1, rec.vis, by = "VisitId")
kn2 = fn1[fn1$Park_Name == nm1,]
kn2$StationID = as.numeric(kn2$StationID)
kn3 <- as.data.frame(table(kn2$StationID,kn2$VisCluster))
names(kn3)[3] <- "Frequency"
p <- ggplot(kn3, aes(Var1,Var2))
p+geom_point(aes(alpha=Frequency,size=Frequency),colour = "Red")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Stations",y="Pattern (Cluster)") +
  ggtitle(paste0("Bubble plot of visit clusters frequency in ",nm1))+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 8,angle = 90,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave(paste0(nm1,"_plot of clusters and statitions_In.jpg"),width = 8, height = 8)
###PART 5
##Finalcodefile.R must run so that 6 outputs can be use in the following sections
#######################################
########Social Netwrok Analysis########
#######################################
##Before code
#Import data
gcode1 <- read.csv(paste0(path1,"/Before net data code.csv"),header = T)
#Prepare adjacency matrix
gcode1[,1] <- names(gcode1)[-1]
row.names(gcode1) <- gcode1[,1]
gcode1[,1] <- NULL
#Prepare the graph using igraph Package
g1 <- graph.adjacency(as.matrix(gcode1),mode="directed",weighted=TRUE)
#Compute the betweeenness using igraph
cen1 <- igraph::betweenness(g1, v = V(g1), directed = TRUE)
#Top 5 betweenness
gsub(pattern = "X","",names(sort(cen1,decreasing = T)[1:10]))
#Clustering the graph using infomap techniques
set.seed(12345)
cm1 <- cluster_infomap(g1)
#Members of each community
table(cm1$membership)
grclus <- as.data.frame(cbind(cm1$names,cm1$membership))
write.csv(grclus,"Before_net_data_code_clus.csv")
#Find the cluster of most used algorithm
obj1 <- grclus[which(grclus$V1 %in% names(sort(cen1,decreasing = T)[1:10])),]
#Filter the nodes which are in the specific cluster
ver1 <- as.character(grclus[which(grclus$V2 %in% c(7)),1])
#Filter garph based on this nodes
g2 <- gcode1[which(rownames(gcode1) %in% ver1),which(colnames(gcode1) %in% ver1)]
#Create a netwrok object using sna package
net1 <- network(x = g2,vertex.attrnames = row.names(gcode1),directed = T)
set.seed(12345)
png(filename="Before_Freq_code.png")
gplot(net1,gmode="digraph",displaylabels = T,object.scale=0.009,label.cex = 1.1,edge.col = "Blue")
dev.off()
###PART 6
#################
##Before time
#Import data
gcode2 <- read.csv(paste0(path1,"/Before net data time.csv"),header = T)
#Prepare adjacency matrix
gcode2[,1] <- names(gcode2)[-1]
row.names(gcode2) <- gcode2[,1]
gcode2[,1] <- NULL
gcode3 <- gcode2/gcode1
gcode3[is.na(gcode3)] <- 0
#Prepare the graph using igraph Package
g1 <- graph.adjacency(as.matrix(gcode3),mode="directed",weighted=TRUE)
#Compute the betweeenness using igraph
cen1 <- igraph::betweenness(g1, v = V(g1), directed = TRUE)
#Top 5 betweenness
gsub(pattern = "X","",names(sort(cen1,decreasing = T)[1:10]))
#"3130"  "13902" "1001"  "13"    "14"    "8"     "13900" "5122"  "59"    "10105"
#Clustering the graph using infomap techniques
set.seed(12345)
cm1 <- cluster_infomap(g1)
table(cm1$membership)
#Save the data
grclus <- as.data.frame(cbind(cm1$names,cm1$membership))
write.csv(grclus,"Before_net_data_time_clus.csv")
#Find the cluster of most used algorithm
obj1 <- grclus[which(grclus$V1 %in% names(sort(cen1,decreasing = T)[1:10])),]
#Filter the nodes which are in the specific cluster
ver1 <- as.character(grclus[which(grclus$V2 %in% c(2)),1])
#Filter garph based on this nodes
g2 <- gcode1[which(rownames(gcode1) %in% ver1),which(colnames(gcode1) %in% ver1)]
#Create a netwrok object using sna package
net1 <- network(x = g2,vertex.attrnames = row.names(gcode1),directed = T)
set.seed(12345)
png(filename="Before_Freq_time.png")
gplot(net1,gmode="digraph",displaylabels = T,object.scale=0.009,label.cex = 1.1,edge.col = "Blue")
dev.off()
###PART 7
###########
##After code
#Import data
gcode1 <- read.csv(paste0(path1,"/After net data code.csv"),header = T)
#Prepare adjacency matrix
gcode1[,1] <- names(gcode1)[-1]
row.names(gcode1) <- gcode1[,1]
gcode1[,1] <- NULL
#Prepare the graph using igraph Package
g1 <- graph.adjacency(as.matrix(gcode1),mode="directed",weighted=TRUE)
#Compute the betweeenness using igraph
cen1 <- igraph::betweenness(g1, v = V(g1), directed = TRUE)
#Top 5 betweenness
gsub(pattern = "X","",names(sort(cen1,decreasing = T)[1:10]))
#"3130"  "13902" "13"    "14"    "1001
#Clustering the graph using infomap techniques
set.seed(12345)
cm1 <- cluster_infomap(g1)
table(cm1$membership)
grclus <- as.data.frame(cbind(cm1$names,cm1$membership))
write.csv(grclus,"After_net_data_code_clus.csv")
#Find the cluster of most used algorithm
obj1 <- grclus[which(grclus$V1 %in% names(sort(cen1,decreasing = T)[1:5])),]
#Filter the nodes which are in the specific cluster
ver1 <- as.character(grclus[which(grclus$V2 %in% c(3)),1])
#Filter garph based on this nodes
g2 <- gcode1[which(rownames(gcode1) %in% ver1),which(colnames(gcode1) %in% ver1)]
#Create a netwrok object using sna package
net1 <- network(x = g2,vertex.attrnames = row.names(gcode1),directed = T)
set.seed(12345)
png(filename="After_Freq_code1.png")
gplot(net1,gmode="digraph",displaylabels = T,object.scale=0.009,label.cex = 1.1,edge.col = "Blue")
dev.off()
###PART 8
###########
##After time
#Import data
gcode2 <- read.csv(paste0(path1,"/After net data time.csv"),header = T)
#Prepare adjacency matrix
gcode2[,1] <- names(gcode2)[-1]
row.names(gcode2) <- gcode2[,1]
gcode2[,1] <- NULL
gcode3 <- gcode2/gcode1
gcode3[is.na(gcode3)] <- 0
#Prepare the graph using igraph Package
g1 <- graph.adjacency(as.matrix(gcode3),mode="directed",weighted=TRUE)
#Compute the betweeenness using igraph
cen1 <- igraph::betweenness(g1, v = V(g1), directed = TRUE)
#Top 5 betweenness
gsub(pattern = "X","",names(sort(cen1,decreasing = T)[1:10]))
#"3130"  "13902" "13"    "14"    "1001
#Clustering the graph using infomap techniques
set.seed(12345)
cm1 <- cluster_infomap(g1)
table(cm1$membership)
grclus <- as.data.frame(cbind(cm1$names,cm1$membership))
write.csv(grclus,"After_net_data_time_clus.csv")
#Find the cluster of most used algorithm
obj1 <- grclus[which(grclus$V1 %in% names(sort(cen1,decreasing = T)[1:5])),]
#Filter the nodes which are in the specific cluster
ver1 <- as.character(grclus[which(grclus$V2 %in% c(2)),1])
#Filter garph based on this nodes
g2 <- gcode1[which(rownames(gcode1) %in% ver1),which(colnames(gcode1) %in% ver1)]
#Create a netwrok object using sna package
net1 <- network(x = g2,vertex.attrnames = row.names(gcode1),directed = T)
set.seed(1234)
png(filename="After_Freq_time1.png")
gplot(net1,gmode="digraph",displaylabels = T,object.scale=0.009,label.cex = 1.1,edge.col = "Blue")
dev.off()
###PART 9
###########
##In code
#Import data
gcode1 <- read.csv(paste0(path1,"/In net data code.csv"),header = T)
#Prepare adjacency matrix
gcode1[,1] <- names(gcode1)[-1]
row.names(gcode1) <- gcode1[,1]
gcode1[,1] <- NULL
#Prepare the graph using igraph Package
g1 <- graph.adjacency(as.matrix(gcode1),mode="directed",weighted=TRUE)
#Compute the betweeenness using igraph
cen1 <- igraph::betweenness(g1, v = V(g1), directed = TRUE)
#Top 5 betweenness
gsub(pattern = "X","",names(sort(cen1,decreasing = T)[1:10]))
#"1020"  "1001"  "1005"  "13902" "1023" 
#Clustering the graph using infomap techniques
set.seed(12345)
cm1 <- cluster_infomap(g1)
table(cm1$membership)
grclus <- as.data.frame(cbind(cm1$names,cm1$membership))
write.csv(grclus,"In_net_data_code_clus.csv")
#Find the cluster of most used algorithm
obj1 <- grclus[which(grclus$V1 %in% names(sort(cen1,decreasing = T)[1:5])),]
#Filter the nodes which are in the specific cluster
ver1 <- as.character(grclus[which(grclus$V2 %in% c(2)),1])
#Filter garph based on this nodes
g2 <- gcode1[which(rownames(gcode1) %in% ver1),which(colnames(gcode1) %in% ver1)]
#Create a netwrok object using sna package
net1 <- network(x = g2,vertex.attrnames = row.names(gcode1),directed = T)
set.seed(1234)
png(filename="In_Freq_code.png")
gplot(net1,gmode="digraph",displaylabels = T,object.scale=0.009,label.cex = 1.1,edge.col = "Blue")
dev.off()
###PART 10
###########
##In time
#Import data
gcode2 <- read.csv(paste0(path1,"/In net data time.csv"),header = T)
gcode2[,1] <- names(gcode2)[-1]
row.names(gcode2) <- gcode2[,1]
gcode2[,1] <- NULL
gcode3 <- gcode2/gcode1
gcode3[is.na(gcode3)] <- 0
#Prepare the graph using igraph Package
g1 <- graph.adjacency(as.matrix(gcode3),mode="directed",weighted=TRUE)
#Compute the betweeenness using igraph
cen1 <- igraph::betweenness(g1, v = V(g1), directed = TRUE)
#Top 5 betweenness
gsub(pattern = "X","",names(sort(cen1,decreasing = T)[1:10]))
#"8"     "13902" "9"     "1023"  "1001"  "1020"  "2"     "7"     "18"    "14001"
#Clustering the graph using infomap techniques
set.seed(12345)
cm1 <- cluster_infomap(g1)
table(cm1$membership)
grclus <- as.data.frame(cbind(cm1$names,cm1$membership))
write.csv(grclus,"In_net_data_time_clus.csv")
#Find the cluster of most used algorithm
obj1 <- grclus[which(grclus$V1 %in% names(sort(cen1,decreasing = T)[1:5])),]
#Filter the nodes which are in the specific cluster
ver1 <- as.character(grclus[which(grclus$V2 %in% c(3)),1])
#Filter garph based on this nodes
g2 <- gcode1[which(rownames(gcode1) %in% ver1),which(colnames(gcode1) %in% ver1)]
#Create a netwrok object using sna package
net1 <- network(x = g2,vertex.attrnames = row.names(gcode1),directed = T)
set.seed(1234)
png(filename="In_Freq_time.png")
gplot(net1,gmode="digraph",displaylabels = T,object.scale=0.009,label.cex = 1.1,edge.col = "Blue")
dev.off()
###PART 11
#####################################################
##Park clustering based on visitIds for visualization
#This process will give the clusters and the order of parks based on the cluster
#First you have run PART 1 until line 160 which is mentioned in the comment
fn1 <- unique(dt1[,which(names(dt1) %in% vars1)])
fn1 <- merge(fn1,rec.vis,by="VisitId")
tb1 <- as.data.frame(table(fn1$Park_Name,fn1$VisCluster))
tb2 <- as.data.frame.matrix(table(fn1$Park_Name,fn1$VisCluster))
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,9)
mat1 <- as.data.frame(matrix(tb1$row1,nrow =37 ,ncol=9))
rownames(mat1) <- sort(unique(dt1$Park_Name))
colnames(mat1) <- 1:9
pcadt.sd <- PCfun(mat1,pcnum = 5)
set.seed(1357)
km1 <- kmeans(pcadt.sd$dt,centers = 9)
mrg1 <- cbind(mat1,km1$cluster)
mrg1$Var1 <- rownames(mrg1)
mrg1 <- mrg1[order(mrg1$`km1$cluster`),]
or1 <- mrg1$Var1
###Plot of park clusters
tb1 <- as.data.frame(table(fn1$Park_Name,fn1$VisCluster))
tb2 <- as.data.frame.matrix(table(fn1$Park_Name,fn1$VisCluster))
tb1$Var1 <- factor(tb1$Var1 , levels = or1)
rs1 <- rowSums(as.matrix(tb2))
tb1$row1 <- tb1$Freq/rep(rs1,9)
rs2 <- colSums(as.matrix(tb2))
rep3=c()
for (i in 9:1){
  obj1=rep(rs2[i],37)
  rep3=  append(obj1,rep3)
}
tb1$col1 <- tb1$Freq/rep3
tb1$newfreq <- discretize(X = tb1$Freq, disc="equalfreq",nbins = 5)
tb1 <- merge(tb1,mrg1,by="Var1")
tb1$`km1$cluster` <- as.factor(tb1$`km1$cluster`)
p <- ggplot(tb1, aes(Var2,Var1))
p+geom_point(aes(alpha=row1,size=row1,colour=`km1$cluster`) )+scale_size_continuous(range = c(0.5,8))+
  labs(x="Clusters",y="Parks",size="Frequency",colour="Cluster",alpha="Frequency") +
  ggtitle("Bubble plot of visit patterns(clusters) - Normalized by Park")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("Clustred Parks based on VisitId.jpg",width = 8, height = 8)
###PART 12
####################
##Clustering Parks##
mlt1 <-  melt(dt1,id.vars =  names(dt1)[c(2:21)])
dt2 <- dcast(mlt1, Park_Name~ value, value.var = c("Park_Name"),fun.aggregate = length)
dt3.std <- dt2[,-1]
for (i in 1:ncol(dt3.std)){
  dt3.std[,i] <- infotheo::discretize(X = dt3.std[,i], disc="equalfreq",nbins = 5)
}
#Run the K-means - Records are VisitIds
d1 <- dist(dt3.std)
h1 <- hclust(d1,method = "ward.D")
plot(h1)
ct1 <- cutree(h1,k=7)
table(ct1)
parkclus <- as.data.frame(cbind(as.character(dt2$Park_Name),ct1))
names(parkclus) <- c("Park","ParkCluster")
###PART 13
#####################################
#####Cluster the Codes###############
mlt1 = melt(dt1,id.vars =  names(dt1)[c(2:21)])
dt2 <- dcast(mlt1, VisitId  ~ value, value.var = c("VisitId"),fun.aggregate = length)
dt3.std <- dt2[,-1]
for (i in 1:ncol(dt3.std)){
  set.seed(12345)
  dt3.std[,i] <- infotheo::discretize(X = dt3.std[,i], disc="globalequalwidth",nbins = 5)
}
pcadt1.zo <- PCfun(t(dt3.std),pcnum = 7)$dt
##First clustering
set.seed(12345)
d <- dist(pcadt1.zo, method = "euclidean")
hclus1 <- hclust(d, method="ward.D") 
plot(hclus1)
groups <- cutree(hclus1, k=5) 
table(groups)
##Second clustering on cluster 2 with 432 members
d <- dist(pcadt1.zo[which(groups==1),], method = "euclidean")
hclus1 <- hclust(d, method="ward.D")
plot(hclus1)
groups2 <- cutree(hclus1, k=5) 
table(groups2)
##Data creation for code clustering
pr1 <- as.data.frame(cbind(names(groups[which(groups!=1)]),unname(groups[which(groups!=1)])))
pr2 <- as.data.frame(cbind(names(groups2),unname(groups2)))
pr2$V2 <- recode(pr2$V2,"2=6;3=7;4=8;5=9")
rec.code <- as.data.frame(rbind(pr1,pr2))
names(rec.code) <- c("Code","CodeCluster")
table(rec.code$CodeCluster)
write.csv(rec.code,"reccode.csv",row.names = FALSE)
##Whether manual stop codes exist on
dtvcod <- merge(dt1,rec.code,by="Code")
rec.code[which(rec.code$Code %in% dt5[dt5$IsManualStop.=="TRUE",][,1]),]
table(rec.code$CodeCluster)
#Table of code with rspect to clusters
rec2 <- merge(rec.code,dt5,by="Code")
tb2 <- as.data.frame.matrix(table(rec2$StopUrgency,rec2$CodeCluster))
#Graph of code frequency with rspect to clusters
tb2 <- as.data.frame.matrix(table(dtvcod$StopUrgency,dtvcod$CodeCluster))
vec1 <- as.vector(as.matrix(tb2))
Var1 <- rep(c(0:6),9)
Var2 <- c(rep(2,7),rep(3,7),rep(4,7),rep(5,7),rep(1,7),rep(6,7),rep(7,7),rep(8,7),rep(9,7))
tb1 <- as.data.frame(cbind(Var1,Var2,vec1))
tb1$col1 <- as.vector(as.matrix(as.data.frame(lapply(tb2, norm1))))
tb1$row1 <- as.vector(as.matrix(t(as.data.frame(lapply(as.data.frame(t(tb2)), norm1)))))
tb1$Var1 <- as.factor(tb1$Var1)
tb1$Var2 <- as.factor(tb1$Var2)
p <- ggplot(tb1, aes(Var2,Var1))+geom_point(aes(alpha=row1,size=row1),colour = "darkGreen")+scale_size_continuous(range = c(0.5,8))+
  labs(x="Code Cluster",y="Urgency",size="Frequency",alpha="Frequency") +
  ggtitle("Bubble plot of code clusters and urgency clusters - Normalized by Urgency")+
  theme(axis.text.x = element_text(hjust=1,vjust=0.5,size = 12,face = "bold",color = "Black"),
        plot.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        axis.title.x =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        axis.title.y =element_text(hjust=0.5,size = 16,face = "bold",color = "Black"),
        legend.title = element_text(hjust=0.5,size = 14,face = "bold",color = "Black"),
        legend.text = element_text(hjust=0.5,size = 10,face = "bold",color = "Black"))
ggsave("Urgency_Code cluster.jpg",width = 10, height = 10)