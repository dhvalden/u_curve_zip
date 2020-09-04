
install.packages('lme4')

library(lme4)
library(ggplot2)
library(mice)


dat <- read.csv('data/fulldataset.csv', row.names = 1)
names(dat)
str(dat)

dat$polori[dat$polori == 77] <- NA
dat$polori[dat$polori == 66] <- NA
## adding gini and gdp in missing cases for 2016
## with closest year or other sources

unique(dat$Sample)

ginitable <- aggregate(gini2016 ~ Sample, unique, data = dat)
ginitable <- rbind(ginitable, list("Australia", 34.4))
ginitable <- rbind(ginitable, list("Canada", 33.8))
ginitable <- rbind(ginitable, list("Chile", 44.4))
ginitable <- rbind(ginitable, list("Czech", 25.4))
ginitable <- rbind(ginitable, list("Russia", 36.8))
str(ginitable)


gdptable <- aggregate(gdp2016 ~ Sample, mean, data = dat)
gdptable <- rbind(gdptable, list("Czech", 18463.3866))
gdptable <- rbind(gdptable, list("Russia", 8704.8984))
str(gdptable)


dat$gini2016 <- NULL
dat$gdp2016 <- NULL

dat <- merge(dat, ginitable, by = 'Sample')
dat <- merge(dat, gdptable, by = 'Sample')
summary(dat)

## trimming down the dataset

dat <- dat[c('Sample',
             'gen',
             'age',
##             'polori',
             'grpid',
             'gai',
             'gbgr',
             'wilac',
             'colac',
             'ins_stigma',
             'perc_stigma',
             'gini2016',
             'gdp2016')]

str(dat)

## temporal multiple imputation

imi <- mice(dat, m = 5, maxit  = 10, seed = 1234)

dat <- complete(imi)
summary(dat)
str(dat)

## auxiliary function

get_mean_centered <- function(y, x, data){
    if (!is.character(y) | !is.character(x)){
        stop('Variables must be entered as text')
    }
    if (!is.data.frame(data)){
        stop('data must be of class = data.frame')
    }
    
    gm_table <- aggregate(data[y], by = data[x], FUN = 'mean',
                          na.rm = TRUE)
    colnames(gm_table) <- c(x, paste0(y, '_gm'))
    tempdata <- merge(data, gm_table, by = x)
    tempdata[paste0(y, '_gmc')] <- tempdata[y] - tempdata[paste0(y, '_gm')]
    tempdata[paste0(y, '_mc')] <- tempdata[y] - mean(tempdata[, y],
                                                     na.rm = TRUE)
    return(tempdata)
}

## creating centered variables
dat <- get_mean_centered('perc_stigma', 'Sample', dat)
dat <- get_mean_centered('age', 'Sample', dat)
## dat <- get_mean_centered('polori', 'Sample', dat)
dat <- get_mean_centered('grpid', 'Sample', dat)

str(dat)
names(dat)

## null models
summary(lmer(wilac ~ 1 + (1 | Sample), data=dat))
summary(lmer(colac ~ 1 + (1 | Sample), data=dat))

## individual level, random intercepts

fit <- lmer(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                ## control variables
                poly(perc_stigma_gmc, 2, raw=FALSE) +
                as.factor(gen) +
                age_gmc +
                grpid_gmc +
                scale(gini2016) +
                scale(gdp2016) +
                (1 | Sample),
            data = dat)

summary(fit)

ggplot(data=dat,
       aes(x = ins_stigma, 
           y = wilac, 
           col = as.factor(Sample)))+
    geom_point(size     = .7,
               alpha    = .8, 
               position = "jitter")+
    stat_smooth(method = lm,
                formula = y ~ x + I(x^2),
                se     = FALSE,
                size   = 1,
                alpha  = .8)
