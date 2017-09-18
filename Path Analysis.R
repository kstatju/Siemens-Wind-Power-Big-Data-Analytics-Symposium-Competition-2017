


#########################################################################

# Multinomial Logistic Regression Model


# Reference http://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/

ml = data
ml$EventWarningStop = relevel(ml$EventWarningStop, ref = "Stop")
fit <- multinom(EventWarningStop ~ Park_Name + Factor + VisitType + ManualStop.during.Visit + IsManualStop. + StopUrgency, data = ml)
ztest = summary(fit)$coefficients/summary(fit)$standard.errors
p <- (1 - pnorm(abs(ztest), 0, 1)) * 2
expbeta = exp(coef(fit))


