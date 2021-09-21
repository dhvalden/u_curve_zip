library(psych)

df <- read.csv2('data/zipdata2.csv')

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

## correcting gender minority self-classification

gen_min_gen <- clean_df$gen >= 3
gen_min_trans <- clean_df$trans < 4

temp_gen_min_df <- data.frame(gen_min_gen, gen_min_trans) 
gen_min_final <- as.numeric(rowSums(temp_gen_min_df, na.rm = TRUE) > 0)

clean_df$gen_min <- gen_min_final

table(clean_df$gen_min)

new_gen <- NA
new_gen[clean_df$gen_min == 1] <- 3
new_gen[clean_df$gen == 1] <- 1
new_gen[clean_df$gen == 2] <- 2

table(clean_df$gen_min, new_gen)
table(new_gen)

clean_df$new_gen <- new_gen

## sanity checks
dim(clean_df)
table(clean_df$Sample)
length(unique(clean_df$Sample))
str(clean_df)

## Saving as cleanzip.csv

write.csv(clean_df, 'data/cleanzip.csv')
