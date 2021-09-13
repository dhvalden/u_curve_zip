
## Loading the raw data downloaded from world bank
giniraw <- read.csv('data/ginidata.csv', skip = 3)
str(giniraw)

gdpraw <- read.csv('data/gdpdata.csv', skip = 3)
str(gdpraw)

## selecting relavnt columns

ginidf <- giniraw[, c('Country.Name', 'X2016')]
colnames(ginidf) <- c('Sample', 'gini2016')

gdpdf <- gdpraw[, c('Country.Name', 'X2016')]
colnames(gdpdf) <- c('Sample', 'gdp2016')

## merging with gai and gdgr

gai_gbgr_df <- read.csv('data/gai_gbgr_data.csv')
colnames(gai_gbgr_df) <- c('Sample', 'gai', 'gbgr')

str(gai_gbgr_df)

countryleveldf <- merge(gai_gbgr_df, ginidf, by = 'Sample', all.x = TRUE)
countryleveldf <- merge(countryleveldf, gdpdf, by = 'Sample', all.x = TRUE)

## merging with individual level clean data

individualdf <- read.csv('data/cleanzip.csv')
str(individualdf)

finaldf <- merge(individualdf, countryleveldf, by = 'Sample')
str(finaldf)

## adding gini and gdp in missing cases for 2016
## with closest year or other sources

ginitable <- aggregate(gini2016 ~ Sample, unique, data = finaldf)
ginitable <- rbind(ginitable, list("Australia", 34.4))
ginitable <- rbind(ginitable, list("Chile", 44.4))
ginitable <- rbind(ginitable, list("Czech", 25.4))
ginitable <- rbind(ginitable, list("Russia", 36.8))
str(ginitable)


gdptable <- aggregate(gdp2016 ~ Sample, mean, data = finaldf)
gdptable <- rbind(gdptable, list("Czech", 18463.3866))
gdptable <- rbind(gdptable, list("Russia", 8704.8984))
str(gdptable)

## Merging again
finaldf$gini2016 <- NULL
finaldf$gdp2016 <- NULL
finaldf <- merge(finaldf, ginitable, by = 'Sample')
finaldf <- merge(finaldf, gdptable, by = 'Sample')

## Sanity checks
summary(finaldf)
str(finaldf)

write.csv(finaldf, 'data/fulldataset.csv')
