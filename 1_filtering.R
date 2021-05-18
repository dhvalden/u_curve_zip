library(psych)

df <- read.csv2('data/zipdata.csv')

str(df)

## removing attention check fails

fltr <- df$att1 == 2 & df$att2 == 6
fltr[is.na(fltr)] <- FALSE
clean_df <- df[fltr, ]
clean_df$X <- NULL
str(clean_df)

## removing samples with less than 15

sample_size <- as.data.frame(table(clean_df$Sample))
colnames(sample_size) <- c('Sample', 'sample_size')

clean_df <- merge(clean_df, sample_size, by='Sample')
clean_df <- clean_df[clean_df$sample_size >= 15, ]

## sanity checks
dim(clean_df)
table(clean_df$Sample)
length(unique(clean_df$Sample))

## Rowwise missing value counts

rowmissing <- rowSums(is.na(clean_df))/20

table(rowmissing)

## Saving as cleanzip.csv

write.csv(clean_df, 'data/cleanzip.csv')
