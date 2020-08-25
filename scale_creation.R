install.packages('lavaan')
install.packages('psych')

library(lavaan)
library(psych)

rawdata <- read.csv('data/fulldataset.csv')

## CFA fitting willingness and participation in collective actions

## Willingness

wilac <- ' willingness  =~ wilac1 + wilac2 + wilac3 +
                           wilac4 +  wilac6 '
fitw <- cfa(wilac, data=rawdata, estimator="MLR",  missing="ML")
summary(fitw, fit.measures=TRUE, modindices=TRUE)

## Participation

colac <- ' participation  =~ colac2 + colac3 +
                             colac4 + colac6 '
fitc <- cfa(colac, data=rawdata, estimator="MLR",  missing="ML")
summary(fitc, fit.measures=TRUE, modindices=TRUE)

## Cronbach Alpha

## wilac
alpha(rawdata[14:19])

## colac
alpha(rawdata[20:25])
