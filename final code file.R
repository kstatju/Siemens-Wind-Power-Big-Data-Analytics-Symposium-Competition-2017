require(stringr)
require(arules)
require(reshape2)
require(ggplot2)
require(lavaan)
require(semPlot)
require(mnlogit)
require(mlogit)
require(AER)
require(nnet)

###############################################
###  CHANGE THE FILE NAME IN datafilename1 AND datafilename2

datafilename1 = "All Sites Together encoded.csv"
datafilename2 = "Codes and Event Warning Stop classification.csv"

### FOLLOWING FILES ARE DATA FILE THAT CONTAINS CLUSTER OF CODE AND CLUSTER OF VISIT ID  AND
### IT CAN BE OBTAINED FORM THE OTHER R-CODE FILE NAMED "CLUSTER.R"
datafilename3 = "reccode.csv"
datafilename4 = "recvis.csv"


### CHANGE THE DATA PATH IN "datapath"
datapath = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Data/Data/"

### CHANGE THE PATH FOR THE R CODE FILE "WIND COMPETITION FUNCTION.R"

funcpath = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/R Code/"

### CHANGE THE PATH FOR SAVING NEW OUTPUT FOR DATA

outdatapath = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/table/"

### CHANGE THE PATH FOR SVING PLOTS
outplotpath = "C:/Users/ka746940/Desktop/UCF/Competition/Wind competition/Output/Plot/"


#### HOPEFULLY YOU DONT NEED TO CHANGE ANYTHIN BELOW. 
#### MIGHT NEED TO INSTALL SOME R PACKAGE IF THERE IS ANY ERROR DURING THE EXECUTION


source(paste0(funcpath, "Wind Competition function.R", collapse = ""))

data1 = read.csv(paste0(datapath,datafilename1, collapse = ""))
data1['visit_start_date_Time'] = as.POSIXct(data1$VisitStartTime,format = '%m/%d/%Y %H:%M')
data1['Time_On'] = as.POSIXct(data1$TimeOn,format = '%m/%d/%Y %H:%M')
data1['Time_Off'] = as.POSIXct(data1$TimeOff,format = '%m/%d/%Y %H:%M')
data1 = subset( data1, select = -c(VisitStartTime, TimeOn, TimeOff ) )
data2 = read.csv(paste0(datapath,datafilename2, collapse = ""))
data = merge(data1, data2, by="Code")


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

unique_Park_name = unique(data$Park_Name)
unique_StationID = unique(data$StationID)
unique_VisitId = unique(data$VisitId)
unique_visit_start_date = unique(data$visit_start_date)
length(unique_VisitId)
length(unique_visit_start_date)
unique_Code = unique(data$Code)
length(unique_Code)

#############################################
## Code for Grouping Error
# 
# aa = data.frame(table(data$Code))
# gCodeList = aa['Var1'][which(aa$Freq<100),]
# 
# data['new_Code'] = data['Code']
# data['new_Code'][which(data$new_Code %in% gCodeList),] = 9999999
# bb = unique(data$new_Code)

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

write.csv(max5MC_prob, file = paste0(outdatapath, "Conditional Probability.csv", collapse = ""))

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
  ggsave(paste0(outplotpath, "Conditional Probability ", unique_Park_name[i], '.png', collapse = ''))
}

byv = 'Park_Name'
unique_Park_name = unlist(unique(data[byv]))
mcmatx = mcmat(data, timevar = "Time_On", resvar = "EWS_StopUr_B_A_concate", visit = 'VisitId', byvar = byv)

#prob_mat = mcmatx/rowSums(mcmatx)
prob_mat = array(apply(mcmatx, 3, tran_prob_mat), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .001,]

write.csv(max5MC_prob, file = paste0(outdatapath,"Transition Probability.csv", collapse = ''))


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
  ggsave(paste0(outplotpath, "Transition Probability ", unique_Park_name[i], '.png', collapse = ''))
}

prob_mat = array(apply(mcmatx, 3, tran_prob_mat_lift), dim(mcmatx), dimnames = dimnames(mcmatx))

# find best 5 conditional probability
max5MC_prob = max_MC_erro(prob_mat, nmax = 5)
max5MC_prob = max5MC_prob[max5MC_prob$Freq > .0001,]


write.csv(max5MC_prob, file = paste0(outdatapath, "Conditional Lift.csv", collapse = ""))

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
  ggsave(paste0(outplotpath, "Conditional Lift ", unique_Park_name[i], '.png', collapse = ''))
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

write.csv(max5MC_prob, file = paste0(outdatapath,"Conditional Prob of CODE.csv", collapse = ''))


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
  ggsave(paste0(outplotpath, "Conditional Prob of CODE ", unique_Park_name[i], '.png', collapse = ''))
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

write.csv(max5MC_prob, file = paste0(outdatapath, "Conditional Prob of CODE CLUSTER.csv", collapse = ""))

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
  ggsave(paste(outplotpath, "Conditional Prob of CODE CLUSTER ", unique_Park_name[i], '.png', collapse = ''))
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
  ggsave(paste(outplotpath, "Conditional Prob of CODE CLUSTER Time ", unique_Park_name[i], '.png', collapse = ''))
}
###############################################################################

# Common Pattern in different time 

comm_pattern = patBasedOnMC(data, timevar = "Time_On", resvar = "Code", visit = 'VisitId', byvar = 'Code_B_A', firstn = 5, inout = 1:2)
comm_pattern1 = patBasedOnFreq(data, resvar = "Code", visit = 'VisitId', byvar = 'Code_B_A', firstn = 5)
write.csv(comm_pattern, file = paste0(outdatapath,"Common Pattern based on MC.csv", collapse = ''))
write.csv(comm_pattern1, file = paste0(outdatapath,"Common Pattern based on Frequency.csv", collapse = ''))


##########################################################################
## Code Count and Code time data for Network Analysis

mcmatxNT = mcmat(data, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(data, timevar = "Time_On", resvar = "Code", visit = 'VisitId')

bdata = data[data$Code_B_A == '1Before',]
mcmatxNT = mcmat(bdata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(bdata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
write.csv(mcmatxNT, file = paste0(outdatapath,"Before net data code.csv", collapse = ''))
write.csv(mcmatxTimeNT, file = paste0(outdatapath,"Before net data time.csv", collapse = ''))

idata = data[data$Code_B_A == '2In',]
mcmatxNT = mcmat(idata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(idata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
write.csv(mcmatxNT, file = paste0(outdatapath,"In net data code.csv", collapse = ''))
write.csv(mcmatxTimeNT, file = paste0(outdatapath,"In net data time.csv", collapse = ''))

adata = data[data$Code_B_A == '3After',]
mcmatxNT = mcmat(adata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
mcmatxTimeNT = mcmattime(adata, timevar = "Time_On", resvar = "Code", visit = 'VisitId')
write.csv(mcmatxNT, file = paste0(outdatapath,"After net data code.csv", collapse = ''))
write.csv(mcmatxTimeNT, file = paste0(outdatapath,"After net data time.csv", collapse = ''))

##########################################################

mcmatx11 = setNames(melt(mcmatxNT), c('Source', 'Target', 'Weight'))
mcmatx11 = mcmatx11[mcmatx11$Weight > 30, ]
write.csv(mcmatx11, file = paste0(outdatapath,"Network data.csv", collapse = ''), row.names = F)


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
dev.copy(png,paste0(outplotpath, 'Plot/PathPlot EWS_ManualStop.png', collapse = ""))
dev.off()
sumy <-  capture.output(capture.output(summary(model.fit, standardized = TRUE, fit.measures = TRUE)), 
                        file = paste0(outdatapath, 'pathsum_EWS_ManualStop.txt', collapse = ""))

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
dev.copy(png,paste0(outplotpath, 'PathPlot EWS_Urgency.png', collapse = ""))
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit1, standardized = TRUE, fit.measures = TRUE)), 
                        file = paste0(outdatapath, 'pathsum_EWS_Urgency.txt', collapse = ""))

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
dev.copy(png,paste0(outplotpath, 'PathPlot EEWS_B_A.png', collapse = ""))
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit2, standardized = TRUE, fit.measures = TRUE)), 
                        file = paste0(outdatapath, 'pathsum_EWS_B_A.txt', collapse = ""))

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
dev.copy(png,paste0(outplotpath, 'PathPlot Cluster.png', collapse = ""))
dev.off()

sumy <-  capture.output(capture.output(summary(model.fit3, standardized = TRUE, fit.measures = TRUE)), 
                        file =  paste0(outdatapath, 'pathsum_Cluster.txt', collapse = ""))


