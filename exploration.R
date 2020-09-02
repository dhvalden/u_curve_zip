
library(psych)

rawdata <- read.csv('data/fulldataset.csv')

rawdata$X <- NULL
rawdata$X.1 <- NULL
rawdata$X.2 <- NULL

names(rawdata)

table(rawdata$additional_data)

wilac <- aggregate(wilac ~ Sample, data = rawdata, mean)

wilac <- merge(wilac, unique(rawdata[c('Sample', 'gai')]), by='Sample')

plot(wilac$gai, wilac$wilac)




colac <- aggregate(colac ~ Sample, data = rawdata, mean)

colac <- merge(colac, unique(rawdata[c('Sample', 'gai')]), by='Sample')

plot(colac$gai, colac$colac)
