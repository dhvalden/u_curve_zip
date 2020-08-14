
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

finaldf <- merge(individualdf, countryleveldf, by = 'Sample')

write.csv(finaldf, 'data/fulldataset.csv')
