
install.packages('lme4')
install.packages('merTools')
install.packages('Amelia')

library(lme4)
library(ggplot2)
library(Amelia)
library(merTools)


dat <- read.csv('data/fulldataset.csv', row.names = 1)
names(dat)
str(dat)

dat$polori[dat$polori == 77] <- NA
dat$polori[dat$polori == 66] <- NA

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
    ## drop original variable and meam cemtered to avoid colliniarity
    tempdata[paste0(y, '_mc')] <- NULL
    tempdata[y] <- NULL
    return(tempdata)
}

## creating centered variables
dat <- get_mean_centered('perc_stigma', 'Sample', dat)
dat <- get_mean_centered('age', 'Sample', dat)
## dat <- get_mean_centered('polori', 'Sample', dat)
dat <- get_mean_centered('grpid', 'Sample', dat)

str(dat)
names(dat)

## chacking correlations
round(cor(dat[,-c(1)], use = "pairwise"), 3)

## remove colliniar variables
dat$gai <- NULL
dat$gbgr <- NULL

## multiple imputation step

impute.out <- amelia(dat, noms = c('gen'), idvars = c('Sample'), m = 10)
summary(impute.out)

## null models
summary(lmerModList(wilac ~ 1 + (1 | Sample),
                    data=impute.out$imputations))
summary(lmerModList(colac ~ 1 + (1 | Sample),
                    data=impute.out$imputations))

## individual level, random intercepts model

#Willingness by perceive stigma model
wilac_perc <- lmerModList(wilac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = impute.out$imputations)
summary(wilac_perc)


#Participation by perceive stigma model
colac_perc <- lmerModList(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = impute.out$imputations)
summary(colac_perc)

#Willingness by institutional stigma model
wilac_ins <- lmerModList(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = impute.out$imputations)
summary(wilac_ins)


#Participation by institutional stigma model
colac_ins <- lmerModList(colac ~ poly(ins_stigma, 2, raw=FALSE) +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = impute.out$imputations)
summary(colac_ins)


## graphical explorations

dat_imp <- impute.out$imputations$imp1

fit_explo <- lm(colac ~ +
                    ## control variables
                    as.factor(gen) +
                    age_gmc +
                    grpid_gmc +
                    scale(gini2016) +
                    scale(gdp2016),
                data = dat_imp)

summary(fit_explo)
dat$residuals <- resid(fit_explo)

fit_res <- lmer(residuals ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                (1 | Sample), data = dat_imp)
summary(fit_res)

ggplot(dat, aes(perc_stigma_gmc, residuals)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) 
