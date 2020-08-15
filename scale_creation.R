install.packages('lavaan')
install.packages('psych')

library(lavaan)
library(psych)

rawdata <- read.csv('data/fulldataset.csv')

## CFA fitting willingness and participation in collective actions

## Willingness

wilac <- ' willingness  =~ wilac1 + wilac2 + wilac3 +
                           wilac4 + wilac5 + wilac6 '
fit <- cfa(wilac, data=rawdata, estimator="MLR",  missing="ML")
summary(fit, fit.measures=TRUE, modindices=TRUE)

## Participation

colac <- ' participation  =~ colac1 + colac2 + colac3 +
                             colac4 + colac5 + colac6 '
fit <- cfa(colac, data=rawdata, estimator="MLR",  missing="ML")
summary(fit, fit.measures=TRUE, modindices=TRUE)

## Cronbach Alpha

## wilac
alpha(rawdata[14:19])

## colac
alpha(rawdata[20:25])
