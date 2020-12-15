#####Libraries
library(dplyr)
library(TraMineR)
library(reshape2)
library(arulesSequences)
library(arules)
library(ggplot2)
path1="D:/job/Siemens/data/"
#####################
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
################
data1 <- read.csv(paste0(path1,"/All Sites Together encoded.csv"),header =T)
data1$VisitStartTime = as.POSIXct(data1$VisitStartTime,format = '%m/%d/%Y %H:%M')
data1$TimeOn = as.POSIXct(data1$TimeOn,format = '%m/%d/%Y %H:%M')
data1$TimeOff = as.POSIXct(data1$TimeOff,format = '%m/%d/%Y %H:%M')
data1$visnew <- paste0(data1$VisitId,"_",data1$TimeOn)
data1["Code_B_A"] = apply(data1[c('VisitStartTime', 'TimeOn', 'VisitDurMinutes')], 1, Code_B_A)
############################
#Combining Codes to the cart
#This function combines codes that are within 'ztime' of each other.
#Also 'crit' here should be one of the "Before/After/In"
#dt2 is the main dataset
srfun1 <- function(dt2,crit,ztime){
  dt1 <- dt2[which(dt2$Code_B_A==crit),]
  df <- dt1 %>%
    arrange(Code) %>%
    unique() %>%
    group_by(VisitId,TimeOn) %>%
    summarise(cart=paste(Code,collapse=";")) %>%
    ungroup()
  #Extract the unique visitIds
  un1 <- unique(df$VisitId)
  z=c()
  for (k in 1:length(un1)){
  y=c()
  j=1
  ln1=length(unique(df[df$VisitId==un1[k],]['TimeOn'])[[1]])-1
  if (ln1<1){
    z[k]=list("1")
  }else{
  for(i in 1:ln1){
      if (unique(df[df$VisitId==un1[k],]['TimeOn'])[[1]][i+1] <= unique(df[df$VisitId==un1[k],]['TimeOn'])[[1]][i]+ztime*60){
          y[i]=paste0(j)
        }
    if(unique(df[df$VisitId==un1[k],]['TimeOn'])[[1]][i+1] > unique(df[df$VisitId==un1[k],]['TimeOn'])[[1]][i]+ztime*60){
      j=j+1
      y[i]=paste0(j)
      }
  }
  z[k] <- list(append(1,y))
  print(list(k," out of ",length(un1)))
  }
  }
  df$newz <- unlist(z)
  df$visn1 <- paste0(df$VisitId,"_",df$newz)
  df3 <- df %>%
    arrange(cart) %>%
    unique() %>%
    group_by(visn1) %>%
    summarise(cart1=paste(cart,collapse=";")) %>%
    ungroup()
  df3$visid <- gsub(pattern ="_(.*)" ,replacement ="" ,x = df3$visn1)
  df3$seq1 <- gsub(pattern ="(.*)_" ,replacement ="" ,x = df3$visn1)
  for (i in 1:nrow(df3)){
    df3$cart2[i] <- lapply(list(unique(sort(as.numeric(unlist(strsplit(df3$cart1[i],";")))))),paste, collapse = ";")
  }
  df3$cart2 <- gsub(";",",",df3$cart2)
  df3$items <- paste0("{",df3$cart2,"}")
  df3$SIZE <- count.fields(textConnection(df3$cart2),sep = ",")
  df3$visid <- as.integer(df3$visid)
  df3$seq1 <- as.integer(df3$seq1)
  df3 <- df3[(order(df3$visid,df3$seq1)),]
  tr1 <- as.data.frame(df3$cart2)
  df6 <- as(tr1,"transactions")
  transactionInfo(df6)$sequenceID <- df3$visid
  transactionInfo(df6)$eventID<- df3$seq1
  transactionInfo(df6)$SIZE<- df3$SIZE
  return(df6)
}
#Concatenat codes which are happened within 5 minutes of each other in before dataset
#For 5 mins extract the data
trbefore <- srfun1(data1,"Before",5)
#Convert extracted data to data frame
dttrbefore <- as(trbefore,"data.frame")
inspect(head(trbefore))
#Run CSpade algorithm (support here is 0.01 which can be changed)
fp_before <-cspade(trbefore,parameter = list(support = 0.01))
#Convert extracted sequential rules to data frame
fp_befored <- as(fp_before,"data.frame")

#After dataset for 5 mins
trafter <- srfun1(data1,"After",5)
dttrafter <- as(trafter,"data.frame")
inspect(head(trafter))
fp_after <-cspade(trafter,parameter = list(support = 0.005))
inspect(sort(fp_after, by = "support"))
fp_afterd <- as(fp_after,"data.frame")

#In dataset for 5 mins
trin <- srfun1(data1,"In",5)
dttrin <- as(trin,"data.frame")
inspect(head(trin))
fp_in <-cspade(trin,parameter = list(support = 0.001))
inspect(sort(fp_in, by = "support"))
fp_ind <- as(fp_in,"data.frame")

#Before dataset for 15 mins
trbefore1 <- srfun1(data1,"Before",15)
dttrbefore1 <- as(trbefore1,"data.frame")
inspect(head(trbefore1))
fp_before1 <-cspade(trbefore1,parameter = list(support = 0.005))
inspect(sort(fp_before1, by = "support"))
fp_before1d <- as(fp_before1,"data.frame")

#After dataset for 15 mins
trafter1 <- srfun1(data1,"After",15)
dttrafter1 <- as(trafter1,"data.frame")
inspect(head(trafter1))
fp_after1 <-cspade(trafter1,parameter = list(support = 0.005))
inspect(sort(fp_after1, by = "support"))
fp_after1d <- as(fp_after1,"data.frame")

#In dataset for 15 mins
trin1 <- srfun1(data1,"In",15)
dttrin1 <- as(trin1,"data.frame")
inspect(head(trin1))
fp_in1 <-cspade(trin1,parameter = list(support = 0.005))
inspect(sort(fp_in1, by = "support"))
fp_in1d <- as(fp_in1,"data.frame")

#Before dataset for 30 mins
trbefore2 <- srfun1(data1,"Before",30)
dttrbefore2 <- as(trbefore2,"data.frame")
inspect(head(trbefore2))
fp_before2 <-cspade(trbefore2,parameter = list(support = 0.005))
inspect(sort(fp_before2, by = "support"))
fp_before2d <- as(fp_before2,"data.frame")

#After dataset for 30 mins
trafter2 <- srfun1(data1,"After",30)
dttrafter2 <- as(trafter2,"data.frame")
inspect(head(trafter2))
fp_after2 <-cspade(trafter2,parameter = list(support = 0.005))
inspect(sort(fp_after2, by = "support"))
fp_after2d <- as(fp_after2,"data.frame")

#In dataset for 30 mins
trin2 <- srfun1(data1,"In",30)
dttrin2 <- as(trin2,"data.frame")
inspect(head(trin2))
fp_in2 <-cspade(trin2,parameter = list(support = 0.005))
inspect(sort(fp_in2, by = "support"))
fp_in2d <- as(fp_in2,"data.frame")

#Before dataset for 1 min
trbefore3 <- srfun1(data1,"Before",1)
dttrbefore3 <- as(trbefore3,"data.frame")
inspect(head(trbefore2))
fp_before3 <-cspade(trbefore3,parameter = list(support = 0.05))
inspect(sort(fp_before3, by = "support"))
fp_before3d <- as(fp_before3,"data.frame")

#After dataset for 1 min
trafter3 <- srfun1(data1,"After",1)
dttrafter3 <- as(trafter3,"data.frame")
inspect(head(trafter3))
fp_after3 <-cspade(trafter3,parameter = list(support = 0.005))
inspect(sort(fp_after3, by = "support"))
fp_after3d <- as(fp_after3,"data.frame")

#In dataset for 1 min
trin3 <- srfun1(data1,"In",1)
dttrin3 <- as(trin3,"data.frame")
inspect(head(trin3))
fp_in3 <-cspade(trin3,parameter = list(support = 0.005))
inspect(sort(fp_in3, by = "support"))
fp_in3d <- as(fp_in3,"data.frame")

#Before dataset for less than 1 min
trbefore4 <- srfun1(data1,"Before",0.01)
fp_before4 <-cspade(trbefore4,parameter = list(support = 0.02))
fp_before4d <- as(fp_before4,"data.frame")

#After dataset for less than 1 min
trafter4 <- srfun1(data1,"After",0.01)
fp_after4 <-cspade(trafter4,parameter = list(support = 0.01))
fp_after4d <- as(fp_after4,"data.frame")

#In dataset for less than 1 min
trin4 <- srfun1(data1,"In",0.01)
fp_in4 <-cspade(trin4,parameter = list(support = 0.005))
fp_in4d <- as(fp_in4,"data.frame")

#Compare the results of after and before rules for time span of 15 minutes
af1 <- round(fp_after1d[which(fp_after1d$sequence %in% intersect(fp_after1d$sequence,fp_before1d$sequence)),2],3)
bf1 <- round(fp_before1d[which(fp_before1d$sequence %in% intersect(fp_after1d$sequence,fp_before1d$sequence)),2],3)
rl1 <- intersect(fp_after1d$sequence,fp_before1d$sequence)
#comp1 contains the table with the same rule in before and after with thier support
comp1 <- as.data.frame(cbind(rl1,af1,bf1))
comp1$rl1 <- gsub("df3\\$cart2\\=","",comp1$rl1)
comp1$rl1 <- gsub("<","",comp1$rl1)
comp1$rl1 <- gsub(">","",comp1$rl1)
# Generate data
bar1 <- ggplot(comp1, aes(factor(af1)))
bar1 + geom_bar()
