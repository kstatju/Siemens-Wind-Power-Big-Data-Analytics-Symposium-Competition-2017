
##############  Required library to perform LSD   ###################
##################################################################

library(agricolae)

################### ONE WAY ANOVA ###############################
#################################################################

#### mean_dif is a function to perform ANOVA  ########

mean_dif <- function(Num_var,Cat_var){       ## Num_var = Numerical variable  Cat_var = Categorical variable
model<-lm(Num_var~Cat_var)
Anova_test <- anova(model)

return(Anova_test)

}


### ANOVA : Average visit time duration in several categorical variables

mean_dif(Pattern_data$VisitDurMinutes,Pattern_data$Park_Name)
mean_dif(Pattern_data$VisitDurMinutes,Pattern_data$StationID)
mean_dif(Pattern_data$VisitDurMinutes,Pattern_data$FactorA)
mean_dif(Pattern_data$VisitDurMinutes,Pattern_data$FactorB)
mean_dif(Pattern_data$VisitDurMinutes,Pattern_data$FactorC)
mean_dif(Pattern_data$VisitDurMinutes,Pattern_data$FactorD)
mean_dif(Pattern_data$VisitDurMinutes,Pattern_data$VisitType)


### ANOVA : Average Time difference in several categorical variables
# Time_diff = Time difference is a new variable that is difined in the report.

mean_dif(Pattern_data$Time_Diff,Pattern_data$Park_Name)
mean_dif(Pattern_data$Time_Diff,Pattern_data$StationID)
mean_dif(Pattern_data$Time_Diff,Pattern_data$FactorA)
mean_dif(Pattern_data$Time_Diff,Pattern_data$FactorB)
mean_dif(Pattern_data$Time_Diff,Pattern_data$FactorC)
mean_dif(Pattern_data$Time_Diff,Pattern_data$FactorD)
mean_dif(Pattern_data$Time_Diff,Pattern_data$VisitType)




##### Least Significant Difference (LSD) post hoc test #########
## Least function is created to perform LSD ####


Least<- function(Num_var,Cat_var){
  LSD<-LSD.test(lm(Num_var~Cat_var), "Cat_var", alpha = 0.05)
  return(LSD$groups)
}



#########  LSD for visit duration   ##############

Least(Pattern_data$VisitDurMinutes,Pattern_data$Park_Name)
Least(Pattern_data$VisitDurMinutes,Pattern_data$StationID)
Least(Pattern_data$VisitDurMinutes,Pattern_data$FactorA)
Least(Pattern_data$VisitDurMinutes,Pattern_data$FactorC)
Least(Pattern_data$VisitDurMinutes,Pattern_data$FactorD)
Least(Pattern_data$VisitDurMinutes,Pattern_data$VisitType)


######  LSD for Time difference ########
### Time difference variable has been explained in the report  ######


Least(Pattern_data$Time_Diff,Pattern_data$Park_Name)
Least(Pattern_data$Time_Diff,Pattern_data$StationID)
Least(Pattern_data$Time_Diff,Pattern_data$FactorA)
Least(Pattern_data$Time_Diff,Pattern_data$FactorC)




###### Boxplot ################
### This boxplot is used to verify the ANOVA test resuts of mean visit duration in parks   ################

####  levels option is used to have a desired order as we got from LSD (visit duration in parks) #############
Patterns2$Park_Name <- factor(Patterns2$Park_Name, levels=c("Park036", "Park015", "Park003","Park006","Park014", "Park007","Park025","Park030","Park026","Park035","Park011","Park029","Park004","Park005","Park009","Park022","Park016","Park019","Park020","Park023","Park021","Park010","Park013","Park027","Park034","Park017","Park002","Park033","Park018","Park012","Park008","Park031","Park037","Park032","Park028","Park001","Park024"))


boxplot(Patterns2$VisitDurMinutes~Patterns2$Park_Name, outline=FALSE, las=2, col="gray", ylab="Visit Duration")