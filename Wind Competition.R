require(stringr)
require(arules)
require(reshape2)
require(ggplot2)
datafilename1 = "All Sites Together encoded.csv"
datafilename2 = "Codes and Event Warning Stop classification.csv"
datafilename3 = "reccode.csv"
datafilename4 = "recvis.csv"

datapath = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/"
funcpath = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/R Code/"
source(paste0(funcpath, "Wind Competition function.R", collapse = ""))

data1 = read.csv(paste0(datapath,datafilename1, collapse = ""))
data1['visit_start_date_Time'] = as.POSIXct(data1$VisitStartTime,format = '%m/%d/%Y %H:%M')
data1['Time_On'] = as.POSIXct(data1$TimeOn,format = '%m/%d/%Y %H:%M')
data1['Time_Off'] = as.POSIXct(data1$TimeOff,format = '%m/%d/%Y %H:%M')
data1 = subset( data1, select = -c(VisitStartTime, TimeOn, TimeOff ) )
data2 = read.csv(paste0(datapath,datafilename2, collapse = ""))
data = merge(data1, data2, by="Code")


###################################
data8 = read.csv("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/Code Group.csv")
data = merge(data, data8, by="Code")

location_assets = read.csv("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/Location List with Number of Assets per Location.csv")
names(location_assets) = c('Park_Name', 'Assets')
park_asset = data.frame(table(data[!duplicated(data$StationID), c('Park_Name', 'StationID')]['Park_Name']))
names(park_asset) = c('Park_Name', 'Freq')
location_assets = merge(location_assets, park_asset, by = 'Park_Name')
rm(park_asset)

data3 = read.csv(paste0(datapath,datafilename3, collapse = ""))
data = merge(data, data3, by="Code")


data4 = read.csv(paste0(datapath,datafilename4, collapse = ""))
data = merge(data, data4, by="VisitId")
#################################################################
data$IsManualStop. = NULL
rm(data2,data1, data3, data4)


data['visit_start_date'] = format(data$visit_start_date_Time,"%d/%m/%Y")
data['visit_start_Time'] = format(data$visit_start_date_Time,"%H:%M")
data['Time_On_date'] = format(data$Time_On,"%d/%m/%Y")
data['Time_On_Time'] = format(data$Time_On,"%H:%M")
data['Time_Off_date'] = format(data$Time_Off,"%d/%m/%Y")
data['Time_Off_Time'] = format(data$Time_Off,"%H:%M")


##################################
## Code for creating variable represents Error occurring "Before", "On", "After" Visit

data['Code_B_A'] = apply(data[c('visit_start_date_Time', 'Time_On', 'VisitDurMinutes')], 1, Code_B_A)

#################################


############################################
## concate EventWarningStop, ManualStop and StopUrgency Variables

data['Factor'] = apply(data[c('FactorA', 'FactorB', 'FactorC', 'FactorD')], 1, concatvar)
data['Code_B_A_concate'] = apply(data[c('Code', 'Code_B_A')], 1, concatvar)
data['EWS_Urgency_concate'] = apply(data[c('EventWarningStop', 'StopUrgency')], 1, concatvar)
data['EWS_Urgency_MStop_concate'] = apply(data[c('EventWarningStop', 'StopUrgency', 'ManualStop')], 1, concatvar)
data['EWS_ManualStop_concate_time'] = apply(data[c('EventWarningStop', 'ManualStop')], 1, concatvar)
data['EWS_B_A_concate'] = apply(data[c('EventWarningStop', 'Code_B_A')], 1, concatvar)
data['EWS_StopUr_B_A_concate'] = apply(data[c('Code_B_A', 'EventWarningStop','StopUrgency')], 1, concatvar)
data['EWS_Urg_MStop_B_A__concate'] = apply(data[c('EventWarningStop', 'StopUrgency', 'ManualStop', 'Code_B_A')], 1, concatvar)
data['Factor'] = apply(data[c('FactorA', 'FactorB','FactorC', 'FactorD')], 1, concatvar)
data['CodeCls_B_A_concate'] = apply(data[c('CodeCluster', 'Code_B_A')], 1, concatvar)


#############################################

# Some Unique factors in variables

names(data)
unique_Park_name = unique(data$Park_Name)
unique_StationID = unique(data$StationID)
unique_VisitId = unique(data$VisitId)
unique_visit_start_date = unique(data$visit_start_date)
length(unique_VisitId)
length(unique_visit_start_date)

unique_Code = unique(data$Code)


#############################################
## Code for Grouping Error

aa = data.frame(table(data$Code))
gCodeList = aa['Var1'][which(aa$Freq<100),]

data['new_Code'] = data['Code']
data['new_Code'][which(data$new_Code %in% gCodeList),] = 9999999
bb = unique(data$new_Code)

########################################################
## Code for Markov chain matrix for Variable EWS_StopUr_B_A_concate

#data1 = data[data$Park_Name == 'Park024',]
byv = 'Park_Name'
unique_Park_name = unlist(unique(data[byv]))
mcmatx = mcmat(data, timevar = "Time_On", resvar = "EWS_StopUr_B_A_concate", visit = 'VisitId', byvar = byv)

#prob_mat = mcmatx/rowSums(mcmatx)
prob_mat = array(apply(mcmatx, 3, tran_prob_mat_all), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .001,]

write.csv(max5MC_prob, file = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/Conditional Probability.csv")

for (i in 1:length(unique_Park_name)){
  a = max5MC_prob[max5MC_prob$Var3 == unique_Park_name[i],]
  a = a[order(a$Var1, a$Var2),]
  xlabel = unique(a$Var1)
  xlabel = xlabel[order(as.character(xlabel))]
  
  p <- ggplot(a, aes(factor(Var1), factor(Var2), size = Freq)) + scale_size(range = c(0, 10)) 
  p + geom_point(aes(colour = factor(xrank))) +scale_x_discrete(limits=xlabel)+scale_y_discrete(limits=xlabel) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste("Conditional Probability Plot of CODE Type, Stop Urgency and \n Time of Code for ", unique_Park_name[i],"\n", collapse = ''),
         x = "Error Type and Time", y = "Error Type and Time", 
         color = "Probability Rank\n", size = 'Probability') +
    theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=8), axis.title.y=element_text(size=9),
          plot.title=element_text(size=13, face="bold", color="darkgreen", hjust = 0.5))
  ggsave(paste("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/Conditional Probability ", unique_Park_name[i], '.png', collapse = ''))
}

byv = 'Park_Name'
unique_Park_name = unlist(unique(data[byv]))
mcmatx = mcmat(data, timevar = "Time_On", resvar = "EWS_StopUr_B_A_concate", visit = 'VisitId', byvar = byv)

#prob_mat = mcmatx/rowSums(mcmatx)
prob_mat = array(apply(mcmatx, 3, tran_prob_mat), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .001,]

write.csv(max5MC_prob, file = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/Transition Probability.csv")


for (i in 1:length(unique_Park_name)){
  a = max5MC_prob[max5MC_prob$Var3 == unique_Park_name[i],]
  a = a[order(a$Var1, a$Var2),]
  xlabel = unique(a$Var1)
  xlabel = xlabel[order(as.character(xlabel))]
  
  p <- ggplot(a, aes(factor(Var1), factor(Var2), size = Freq)) + scale_size(range = c(0, 10)) 
  p + geom_point(aes(colour = factor(xrank))) +scale_x_discrete(limits=xlabel)+scale_y_discrete(limits=xlabel) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste("Transition Probability Plot of CODE Type, Stop Urgency and \n Time of Code for ", unique_Park_name[i],"\n", collapse = ''),
         x = "Error Type and Time", y = "Error Type and Time", 
         color = "Probability Rank\n", size = 'Probability') +
    theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=8), axis.title.y=element_text(size=9),
          plot.title=element_text(size=13, face="bold", color="darkgreen", hjust = 0.5))
  ggsave(paste("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/Transition Probability ", unique_Park_name[i], '.png', collapse = ''))
}

prob_mat = array(apply(mcmatx, 3, tran_prob_mat_lift), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .0001,]


write.csv(max5MC_prob, file = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/Conditional Lift.csv")

for (i in 1:length(unique_Park_name)){
  a = max5MC_prob[max5MC_prob$Var3 == unique_Park_name[i],]
  a = a[order(a$Var1, a$Var2),]
  xlabel = unique(a$Var1)
  xlabel = xlabel[order(as.character(xlabel))]
  
  p <- ggplot(a, aes(factor(Var1), factor(Var2), size = Freq)) + scale_size(range = c(0, 10)) 
  p + geom_point(aes(colour = factor(xrank))) +scale_x_discrete(limits=xlabel)+scale_y_discrete(limits=xlabel) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste("Conditional Lift Plot of CODE Type, Stop Urgency and \n Time of Code for ", unique_Park_name[i],"\n", collapse = ''),
         x = "Error Type and Time", y = "Error Type and Time", 
         color = "Lift Rank\n", size = 'Conditional Lift') +
    theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=8), axis.title.y=element_text(size=9),
          plot.title=element_text(size=13, face="bold", color="darkgreen", hjust = 0.5))
  ggsave(paste("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/Conditional Lift ", unique_Park_name[i], '.png', collapse = ''))
}

########################################################
## Code for Markov chain matrix for variable CODE

#data1 = data[data$Park_Name == 'Park024',]
byv = 'Park_Name'
unique_Park_name = unlist(unique(data[byv]))
mcmatx = mcmat(data, timevar = "Time_On", resvar = "Code", visit = 'VisitId', byvar = byv)

#prob_mat = mcmatx/rowSums(mcmatx)
mcmatx1 = apply(mcmatx, 3, function(x) tranc_tran_matrix(x, 10))
prob_mat = array(apply(mcmatx, 3, tranc_tran_matrix), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .25,]

write.csv(max5MC_prob, file = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/Conditional Prob of CODE.csv")


for (i in 1:length(unique_Park_name)){
  a = max5MC_prob[max5MC_prob$Var3 == unique_Park_name[i],]
  a = a[order(a$Var1, a$Var2),]
  xlabel = unique(a$Var1)
  xlabel = xlabel[order(as.character(xlabel))]
  
  p <- ggplot(a, aes(factor(Var1), factor(Var2), size = Freq)) + scale_size(range = c(0, 10)) 
  p + geom_point(aes(colour = factor(xrank))) +scale_x_discrete(limits=xlabel)+scale_y_discrete(limits=xlabel) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste("Conditional Probability Plot of CODE for ", unique_Park_name[i],"\n", collapse = ''), x = "Given Code", y = "Expected Code", 
         color = "Probability Rank\n", size = 'Probability') +
    theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=8), axis.title.y=element_text(size=9),
          plot.title=element_text(size=13, face="bold", color="darkgreen", hjust = 0.5))
  ggsave(paste("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/Conditional Prob of CODE ", unique_Park_name[i], '.png', collapse = ''))
}


########################################################
## Code for Markov chain matrix for variable CodeCluster

#data1 = data[data$Park_Name == 'Park024',]
byv = 'Park_Name'
unique_Park_name = unlist(unique(data[byv]))
mcmatx = mcmat(data, timevar = "Time_On", resvar = "CodeCluster", visit = 'VisitId', byvar = byv)

#prob_mat = mcmatx/rowSums(mcmatx)
prob_mat = array(apply(mcmatx, 3, tran_prob_mat_all), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .008,]

write.csv(max5MC_prob, file = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/Conditional Prob of CODE CLUSTER.csv")

for (i in 1:length(unique_Park_name)){
  a = max5MC_prob[max5MC_prob$Var3 == unique_Park_name[i],]
  a = a[order(a$Var1, a$Var2),]
  xlabel = unique(a$Var1)
  xlabel = xlabel[order(as.character(xlabel))]
  
  p <- ggplot(a, aes(factor(Var1), factor(Var2), size = Freq)) + scale_size(range = c(0, 10)) 
  p + geom_point(aes(colour = factor(xrank))) +scale_x_discrete(limits=xlabel)+scale_y_discrete(limits=xlabel) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste("Conditional Probability Plot of \n CODE CLUSTER for ", unique_Park_name[i],"\n", collapse = ''), x = "Given Code", y = "Expected Code", 
         color = "Probability Rank\n", size = 'Probability') +
    theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=8), axis.title.y=element_text(size=9),
          plot.title=element_text(size=13, face="bold", color="darkgreen", hjust = 0.5))
  ggsave(paste("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/Conditional Prob of CODE CLUSTER ", unique_Park_name[i], '.png', collapse = ''))
}

########################################################

## Code for Markov chain matrix for variable CodeCluster Before After

#data1 = data[data$Park_Name == 'Park024',]
byv = 'Park_Name'
unique_Park_name = unlist(unique(data[byv]))
mcmatx = mcmat(data, timevar = "Time_On", resvar = "CodeCls_B_A_concate", visit = 'VisitId', byvar = byv)

#prob_mat = mcmatx/rowSums(mcmatx)
prob_mat = array(apply(mcmatx, 3, tran_prob_mat_all), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .009,]

for (i in 1:length(unique_Park_name)){
  a = max5MC_prob[max5MC_prob$Var3 == unique_Park_name[i],]
  a = a[order(a$Var1, a$Var2),]
  xlabel = unique(a$Var1)
  xlabel = xlabel[order(as.character(xlabel))]
  
  p <- ggplot(a, aes(factor(Var1), factor(Var2), size = Freq)) + scale_size(range = c(0, 10)) 
  p + geom_point(aes(colour = factor(xrank))) +scale_x_discrete(limits=xlabel)+scale_y_discrete(limits=xlabel) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    labs(title = paste("Conditional Probability Plot of CODE CLUSTER \n with Time for ", unique_Park_name[i],"\n", collapse = ''), x = "Given Code", y = "Expected Code", 
         color = "Probability Rank\n", size = 'Probability') +
    theme(axis.text.x=element_text(size=8), axis.title.x=element_text(size=9),
          axis.text.y=element_text(size=8), axis.title.y=element_text(size=9),
          plot.title=element_text(size=13, face="bold", color="darkgreen", hjust = 0.5))
  ggsave(paste("C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/Conditional Prob of CODE CLUSTER Time ", unique_Park_name[i], '.png', collapse = ''))
}
###############################################################################

# Common Pattern in different time 

comm_pattern = patBasedOnMC(data, timevar = "Time_On", resvar = "Code", visit = 'VisitId', byvar = 'Code_B_A', firstn = 5, inout = 1:2)
comm_pattern1 = patBasedOnFreq(data, resvar = "Code", visit = 'VisitId', byvar = 'Code_B_A', firstn = 5)
write.csv(comm_pattern, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/Common Pattern based on MC.csv')
write.csv(comm_pattern1, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/Common Pattern based on Frequency.csv')


##########################################################################
## Code Count and Code time data for Network Analysis

mcmatxNT = mcmat(data, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(data, timevar = "Time_On", resvar = "Code", visit = 'VisitId')

bdata = data[data$Code_B_A == '1Before',]
mcmatxNT = mcmat(bdata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(bdata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
write.csv(mcmatxNT, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/Before net data code.csv')
write.csv(mcmatxTimeNT, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/Before net data time.csv')

idata = data[data$Code_B_A == '2In',]
mcmatxNT = mcmat(idata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(idata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
write.csv(mcmatxNT, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/In net data code.csv')
write.csv(mcmatxTimeNT, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/In net data time.csv')

adata = data[data$Code_B_A == '3After',]
mcmatxNT = mcmat(adata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(adata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
write.csv(mcmatxNT, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/After net data code.csv')
write.csv(mcmatxTimeNT, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/After net data time.csv')

# 
# paste0(as.vector(names(sort(prob_mat[2,][prob_mat[2,] != 0], decreasing = TRUE)[1:10])), collapse = ', ')
# prob_mat2 = prob_mat %*% prob_mat
# paste0(as.vector(names(sort(prob_mat2[1,][prob_mat2[1,] != 0], decreasing = TRUE)[1:10])), collapse = ', ')
# prob_mat3 = prob_mat2 %*% prob_mat
# paste0(as.vector(names(sort(prob_mat3[1,][prob_mat3[1,] != 0], decreasing = TRUE)[1:10])), collapse = ', ')

##########################################################

mcmatx11 = setNames(melt(mcmatxNT), c('Source', 'Target', 'Weight'))
mcmatx11 = mcmatx11[mcmatx11$Weight > 30, ]
write.csv(mcmatx11, file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/Network data.csv', row.names = F)


###################################################
## Code for creating Pattern Matrix

varlist1 = c('VisitId', 'EWS_ManualStop_concate_time')
patt_mat = pattrn_mat(data[varlist1], varlist = varlist1)
sum(patt_mat)
##################################################################

###################################################################
## Correlation of Pattern

corrmat1 = corrmat
diag(corrmat1) = 0
pattrn = data.frame(apply(corrmat1, 1, cor_based_pattern))
pa = data.frame(table(pattrn))
cor_based_pattern(corrmat[1,])



pattrn_data = pattern5seq(data)

pattern_data1 = pattrn_data
pattern_data1['Pattern_Concate'] = apply(pattern_data1[c('BEr2', 'IEr2', 'IEr3')], 1, concatvar)
aa = data.frame(table(pattern_data1$Pattern_Concate))
aa = data.frame(table(pattrn_data$Pattern))


colSums(is.na(pattrn_data[,10:24]))


#############################
#path Analysis

require(lavaan)
require(semPlot)
require(mnlogit)
require(mlogit)
require(AER)
require(nnet)


# Path Analysis based on EventWarningStop with Manual Stop

varlist1 = c('VisitId', 'EWS_ManualStop_concate_time')
patt_mat = pattrn_mat(data[varlist1], varlist = varlist1)


corrmat = cor(patt_mat)
rname = row.names(corrmat)
rownames(corrmat) = colnames(corrmat) = rname
model = 'Stop.FALSE ~ a1*Event.FALSE + a2 * Warning.FALSE
Warning.FALSE ~ c1 * Event.FALSE  + c2 * Stop.TRUE
Stop.TRUE ~  b1 * Stop.FALSE 
IndEff1 := a2*c1
IndEff2 := b1 * a1
Total1 := a1 + (a2*c1)
'

model.fit = sem(model, sample.cov = corrmat, sample.nobs = nrow(patt_mat))
summary(model.fit, standardized = TRUE, fit.measures = TRUE)
semPaths(model.fit, what = 'mod', whatLabels = 'par', layout = 'tree2', intercepts = FALSE, style = "ram",
         fixedStyle = c('red',1), freeStyle = c('red', 1), 
         title = TRUE, edge.color = 'red', color = 'green', edge.label.cex = .8, nCharNodes = 3, nCharEdges = 3, 
         sizeMan = 6, curvature = 2, sizeInt = 5, nDigits = 3, curvePivot = TRUE)
title('Path Model by  EventWarningStop and ManualStop')
dev.copy(png,'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/PathPlot EWS_ManualStop.png')
dev.off()
sumy <-  capture.output(capture.output(summary(model.fit, standardized = TRUE, fit.measures = TRUE)), 
               file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/pathsum_EWS_ManualStop.txt')

fitMeasures(model.fit)



# Path Analysis based on EventWarningStop with Stop Urgency 

varlist2 = c('VisitId', 'EWS_Urgency_concate')
patt_mat1 = pattrn_mat(data[varlist2], varlist = varlist2)


corrmat1 = cor(patt_mat1)
rname = row.names(corrmat1)
rownames(corrmat1) = colnames(corrmat1) = rname
model1 = 'Stop.5 ~ a1 * Stop.3 +  a2 * Stop.4 + a3 * Warning.0
Warning.0 ~ b1 * Event.0
Stop.1 ~ c1 * Warning.0
Stop.2 ~ d1 * Stop.1
Stop.3 ~ e1 * Stop.2
Stop.4 ~ f1 * Stop.3 + f2 * Warning.0
Event.0 ~ g1 * Stop.5
Warning.0 ~ h1 * Stop.6
Warning.0 ~~ Event.0
IndEff1 := a3 * b1
IndEff2 := c1 * d1
IndEff3 := e1 * a1 * d1
Total := IndEff3 + a3 + a2
IndEff4 := g1 * a3
IndEff5 := g1 * b1
'

prcomp(iris[1:3,1:4])



model.fit1 = sem(model1, sample.cov = corrmat1, sample.nobs = nrow(patt_mat1))
summary(model.fit1, standardized = TRUE, fit.measures = TRUE)
semPaths(model.fit1, what = 'mod', whatLabels = 'par', layout = 'circle', intercepts = FALSE, style = "ram",
         fixedStyle = c('red',1), freeStyle = c('red', 1), 
         title = TRUE, edge.color = 'red', color = 'green', edge.label.cex = .8, nCharNodes = 3, nCharEdges = 3, 
         sizeMan = 6, curvature = 2, sizeInt = 5, nDigits = 3)
title('Path Model by  EventWarningStop and Stop Urgency')
dev.copy(png,'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/PathPlot EWS_Urgency.png')
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit1, standardized = TRUE, fit.measures = TRUE)), 
                        file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/pathsum_EWS_Urgency.txt')

fitMeasures(model.fit1)



# Path Analysis based on EventWarningStop with Time of Visit (Before, Within, After)


varlist3 = c('VisitId', 'EWS_B_A_concate')
patt_mat2 = pattrn_mat(data[varlist3], varlist = varlist3)


corrmat2 = cor(patt_mat2)
rname = row.names(corrmat2)
rownames(corrmat2) = colnames(corrmat2) = rname
model2 = 'Stop.1Before ~ a1 * Event.1Before + a2 * Warning.1Before
Warning.1Before ~ b1 * Event.1Before 
Event.2In ~  c1 * Stop.2In + c2 * Warning.2In
Warning.2In  ~ d1 * Stop.2In + d2 * Warning.1Before
Warning.3After ~ e1 *Stop.3After
Stop.2In ~ f1 * Stop.1Before
IndEff1 := a2 * b1
IndEff2 := d1 * c2
IndEff3 := b1 * d2
IndEff4 := f1 * a2
'

model.fit2 = sem(model2, sample.cov = corrmat2, sample.nobs = nrow(patt_mat2))
summary(model.fit2, standardized = TRUE, fit.measures = TRUE)
semPaths(model.fit2, what = 'mod', whatLabels = 'par', layout = 'spring', intercepts = FALSE, style = "ram",
         fixedStyle = c('red',1), freeStyle = c('red', 1), 
         title = TRUE, edge.color = 'red', color = 'green', edge.label.cex = .8, nCharNodes = 5, nCharEdges = 3, 
         sizeMan = 6, curvature = 3, sizeInt = 5, nDigits = 3)
title('Path Model by  EventWarningStop and Time of Visit \n(Before, Within, After)')
dev.copy(png,'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/PathPlot EEWS_B_A.png')
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit2, standardized = TRUE, fit.measures = TRUE)), 
                        file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/pathsum_EWS_B_A.txt')

fitMeasures(model.fit2)



# Path Analysis based on Code Cluster 


varlist4 = c('VisitId', 'CodeCluster')
patt_mat3 = pattrn_mat(data[varlist4], varlist = varlist4)


corrmat3 = cor(patt_mat3)
rname = row.names(corrmat3)
rownames(corrmat3) = colnames(corrmat3) = paste("Cluster",rname, sep = '')
model3 = 'Cluster1 ~ a1 * Cluster2 + a2 * Cluster3
Cluster2 ~ b1 * Cluster3 
Cluster3 ~  c1 * Cluster2 + c2 * Cluster4 + c3 * Cluster5 + c4 * Cluster8 + c5 * Cluster9
Cluster4  ~ d1 * Cluster8 + d2 * Cluster3
Cluster5 ~ e1 * Cluster9
Cluster6 ~ f1 * Cluster3
Cluster7 ~ g1 * Cluster1
'

model.fit3 = sem(model3, sample.cov = corrmat3, sample.nobs = nrow(patt_mat3))
summary(model.fit3, standardized = TRUE, fit.measures = TRUE)
semPaths(model.fit3, what = 'mod', whatLabels = 'par', layout = 'tree2', intercepts = FALSE, style = "ram",
         fixedStyle = c('red',1), freeStyle = c('red', 1), 
         title = TRUE, edge.color = 'red', color = 'green', edge.label.cex = .8, nCharNodes = 5, nCharEdges = 3, 
         sizeMan = 6, curvature = 5, sizeInt = 5, nDigits = 3)
title('Path Model by  Code Cluster')
dev.copy(png,'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/PathPlot Cluster.png')
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit3, standardized = TRUE, fit.measures = TRUE)), 
                        file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/pathsum_Cluster.txt')




# Path Analysis based on Code Cluster with Time of Visit (Before, Within, After)


varlist5 = c('VisitId', 'CodeCls_B_A_concate')
patt_mat4 = pattrn_mat(data[varlist5], varlist = varlist5)


corrmat4 = cor(patt_mat4)
rname = row.names(corrmat4)
corrmat4 = corrmat4[c('Clus1.1Before', 'Clus2.1Before', 'Clus3.1Before', 'Clus4.1Before', 'Clus5.1Before', 'Clus6.1Before', 'Clus7.1Before', 'Clus8.1Before', 'Clus9.1Before'),]
rownames(corrmat4) = colnames(corrmat4) = paste('Clus', rname, sep = '')

model4 = "
Clus1.1Before ~ Clus2.1Before + Clus3.1Before + Clus7.1Before
Clus2.1Before ~ Clus3.1Before + Clus7.1Before
Clus3.1Before ~ Clus5.1Before + Clus6.1Before + Clus8.1Before + Clus9.1Before
Clus4.1Before ~ Clus8.1Before + Clus3.1Before
Clus5.1Before ~ Clus9.1Before
Clus6.1Before ~ Clus3.1Before
Clus7.1Before ~ Clus1.1Before + Clus2.1Before
"

Clus1.2In ~ Clus2.2In + Clus3.2In + Clus7.2In + Clus8.2In + Clus1.1Before
Clus2.2In ~ Clus3.2In + Clus5.2In + Clus8.2In + Clus2.1Before
Clus3.2In ~ Clus4.2In + Clus5.2In + Clus6.2In + Clus7.2In + Clus9.2In + Clus3.1Before
Clus4.2In ~ Clus8.2In + Clus4.1Before
Clus5.2In ~ Clus9.2In + Clus5.1Before
Clus7.2In ~ Clus7.1Before
Clus9.2In ~ Clus9.1Before
Clus1.3After ~ Clus2.3After + Clus3.3After + Clus5.3After + Clus1.2In + Clus2.2In
Clus2.3After ~ Clus3.3After + Clus4.3After + Clus8.3After + Clus2.2In
Clus3.3After ~ Clus4.3After + Clus5.3After + Clus6.3After + Clus8.3After + Clus9.3After
Clus4.3After ~ Clus8.3After + Clus9.3After
Clus5.3After ~ Clus9.3After + Clus5.2In 
Clus8.3After ~ Clus9.3After
Clus9.3After ~ Clus5.2In
""

model.fit4 = sem(model4, sample.cov = corrmat4, sample.nobs = nrow(patt_mat4))
summary(model.fit4, standardized = TRUE, fit.measures = TRUE)
semPaths(model.fit4, what = 'mod', whatLabels = 'par', layout = 'spring', intercepts = FALSE, style = "ram",
         fixedStyle = c('red',1), freeStyle = c('red', 1), 
         title = TRUE, edge.color = 'red', color = 'green', edge.label.cex = .8, nCharNodes = 5, nCharEdges = 3, 
         sizeMan = 6, curvature = 5, sizeInt = 5, nDigits = 3)
title('Path Model by  Code Cluster Time of Visit (Before, Within, After)')
dev.copy(png,'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/PathPlot Cluster With Time.png')
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit4, standardized = TRUE, fit.measures = TRUE)), 
                        file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/pathsum_Cluster With time.txt')



##################################################################
# Cluster Analysis

Clus_code1 = clust(data, nclust = 6)

#########################################################################################


# Factor Analysis
EWS_StopUr_B_A_concate

varlist1 = c('VisitId', 'Code')
patt_mat = pattrn_mat(data[varlist1], varlist = varlist1)
patt_mat = as.data.frame(patt_mat)
patt_mat[,which(colSums(patt_mat) < 1000)] <- NULL

pca = princomp(patt_mat, cor = T)
summary(pca)
plot(pca)
covv = cor(patt_mat)
covv = cor.smooth(patt_mat)
det(covv)
# Maximum Likelihood Factor Analysis
# entering raw data and extracting 3 factors, 
# with varimax rotation 
fit <- factanal(patt_mat, 5, rotation="varimax")
print(fit, digits=2, cutoff=.3, sort=TRUE)
# plot factor 1 by factor 2 
load <- fit$loadings[,1:2] 
plot(load,type="n") # set up plot 
text(load,labels=names(patt_mat),cex=.7) # add variable names
library(psych)
fit1 <- fa(patt_mat, nfactors=10, fm="mle", rotate="varimax", max.iter = 1000)
print(fit1$loadings, sort = T, cut=.3)
loding = as.data.frame.array(fit1$loadings)
loding[loding < 0.2] = NA
library(nFactors)
ev <- eigen(cor(patt_mat)) # get eigenvalues
ap <- parallel(subject=nrow(patt_mat),var=ncol(patt_mat),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
plotnScree(nS)

#################################

scor = as.data.frame(fit1$scores)
names(scor) = c("F10", "F5", "F1",  "F2",  "F4",  "F9",  "F8",  "F3",  "F6",  "F7")

patt_mat2 = cbind(patt_mat, scor)
require(sem)
patt_mat[,!(names(patt_mat)%in% c('2200', '3118','3123','2', '7', '9'))] = NULL
corrmat2 = cor(patt_mat)
rname = row.names(corrmat2)
rownames(corrmat2) = colnames(corrmat2) = rname
model2 <- specifyModel()
F10 -> 2200, lam27, NA
F10 -> 3118, lam28, NA
F10 -> 3123, lam29, NA
F5 -> 2, lam1, NA
F5 -> 7, lam2, NA
F5 -> 9, lam4, NA
2200 <-> 2200, e1,   NA 
3118 <-> 3118, e2,   NA 
F1 <-> F1, NA,    1 
F2 <-> F2, NA,    1 
F1 <-> F2, F1F2, NA


model.fit2 = sem::sem(model2, S = corrmat2, N = 100)
summary(model.fit2)
# print standardized coefficients (loadings) 
std.coef(model.fit2)


summary(model.fit2, standardized = TRUE, fit.measures = TRUE)
semPaths(model.fit2, what = 'mod', whatLabels = 'par', layout = 'circle', intercepts = FALSE, style = "ram",
         fixedStyle = c('red',1), freeStyle = c('red', 1), 
         title = TRUE, edge.color = 'red', color = 'green', edge.label.cex = .8, nCharNodes = 3, nCharEdges = 3, 
         sizeMan = 6, curvature = 2, sizeInt = 5, nDigits = 3)
title('Path Model by  Code')
dev.copy(png,'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/PathPlot Code.png')
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit1, standardized = TRUE, fit.measures = TRUE)), 
                        file = 'C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/pathsum_EWS_Urgency.txt')

fitMeasures(model.fit1)













library(sem)
mydata.cov <- cov(mydata)
model.mydata <- specifyModel(text = " 
F1 ->  X1, lam1, NA
F1 ->  X2, lam2, NA 
F1 ->  X3, lam3, NA 
F2 ->  X4, lam4, NA 
F2 ->  X5, lam5, NA 
F2 ->  X6, lam6, NA 
X1 <-> X1, e1,   NA 
X2 <-> X2, e2,   NA 
X3 <-> X3, e3,   NA 
X4 <-> X4, e4,   NA 
X5 <-> X5, e5,   NA 
X6 <-> X6, e6,   NA 
F1 <-> F1, NA,    1 
F2 <-> F2, NA,    1 
F1 <-> F2, F1F2, NA
")
mydata.sem <- sem(model.mydata, mydata.cov, nrow(mydata))