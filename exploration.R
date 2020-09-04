
library(psych)
library(ggplot2)

rawdata <- read.csv('data/fulldataset.csv', row.names = 1)

names(rawdata)

## auxiliary dataframes fro plotting

wilac <- aggregate(wilac ~ Sample, data = rawdata, mean)
wilac <- merge(wilac, unique(rawdata[c('Sample', 'ins_stigma')]), by='Sample')

colac <- aggregate(colac ~ Sample, data = rawdata, mean)
colac <- merge(colac, unique(rawdata[c('Sample', 'ins_stigma')]), by='Sample')

## ploting

ggplot(rawdata, aes(perc_stigma, colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) 

ggplot(colac, aes(ins_stigma, colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) 

ggplot(wilac, aes(ins_stigma, wilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)
