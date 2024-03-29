install.packages('lavaan')
install.packages('psych')

library(lavaan)
library(psych)

rawdata <- read.csv('data/fulldataset.csv')

## CFA fitting willingness and participation in collective actions

## Willingness

wilac <- ' willingness  =~ wilac1 + wilac2 + wilac3 '
fitw <- cfa(wilac, data=rawdata, estimator="MLR",  missing="ML")
summary(fitw, fit.measures=TRUE, modindices=TRUE)

## Participation

colac <- ' participation  =~ colac1 + colac2 + colac3 '
fitc <- cfa(colac, data=rawdata, estimator="MLR",  missing="ML")
summary(fitc, fit.measures=TRUE, modindices=TRUE)

## Cronbach Alpha

## wilac
alpha(rawdata[c('wilac1', 'wilac2', 'wilac3')])
## colac
alpha(rawdata[c('colac1', 'colac2', 'colac3')])


## Pecieve institutional stigma

corr.test(rawdata$accep, rawdata$disclo)
rawdata$perc_stigma <- rowMeans(rawdata[c('accep', 'disclo')], na.rm = TRUE)
rawdata$perc_stigma <- 8 - rawdata$perc_stigma # reversing for interpretability

## calculating composite variables

rawdata$wilac <- rowMeans(rawdata[c('wilac1',
                                    'wilac2',
                                    'wilac3')], na.rm = TRUE)

rawdata$colac <- rowMeans(rawdata[c('colac1',
                                    'colac2',
                                    'colac3')], na.rm = TRUE)

## Adding low cost collective acction for apendix comparison

rawdata$wilac_l <- rowMeans(rawdata[c('wilac4',
                                      'wilac5',
                                      'wilac6')], na.rm = TRUE)

rawdata$colac_l <- rowMeans(rawdata[c('colac4',
                                      'colac5',
                                      'colac6')], na.rm = TRUE)

## correlation beetwen measures of intitutional stigma

corr.test(rawdata$gai, rawdata$gbgr)

## based on our pre-resgistration, the correlation is high enough to conbine the 2 indices
## Also reversing and scaling

rawdata$gai_r <- 10 - rawdata$gai # theoretical maximum 
rawdata$gbgr_r <- 100 - rawdata$gbgr
    
gai_s <- scale(rawdata$gai_r)
gai_s
gbgr_s <- scale(rawdata$gbgr_r)
gbgr_s

rawdata$ins_stigma <- rowMeans(data.frame(gai_s, gbgr_s),
                               na.rm = FALSE)

## deleting useles index variables

rawdata$X <- NULL
rawdata$X.1 <- NULL
rawdata$X.2 <- NULL

dim(rawdata)
head(rawdata)
str(rawdata)

## overwirtting old data file

write.csv(rawdata, 'data/fulldataset.csv')
