library(ggplot2)
library(Amelia)
library(sjPlot)
library(psych)
library(gridExtra)
library(mice)
library(miceadds)
library(brms)

dat <- read.csv('data/fulldataset.csv', row.names = 1)
names(dat)
str(dat)

dat$polori[dat$polori == 77] <- NA
dat$polori[dat$polori == 66] <- NA
dat$polori[dat$polori == 8] <- NA

table(dat$polori)

## trimming down the dataset

dat <- dat[c('ID',
             'Sample',
             'gen_min',
             'age',
             'polori',
             'grpid',
             'gai_r',
             'gbgr_r',
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
    return(tempdata)
}

#################
## descriptives
#################

options(browser="/usr/bin/firefox")

tab_df(aggregate(ID ~ Sample, data = dat, FUN = length))
prop.table(table(dat$gen_min))*100

frqdat <- subset(dat, select = -c(ID, Sample, gen_min))
tab_df(describe(frqdat))
tab_corr(frqdat)

## creating centered variables
dat <- get_mean_centered('perc_stigma', 'Sample', dat)
dat <- get_mean_centered('age', 'Sample', dat)
dat <- get_mean_centered('polori', 'Sample', dat)
dat <- get_mean_centered('grpid', 'Sample', dat)
dat$ins_stigma <- as.vector(scale(dat$ins_stigma))
dat$gini2016s <- as.vector(scale(dat$gini2016))
dat$gdp2016s <- as.vector(scale(dat$gdp2016))

str(dat)
names(dat)
summary(dat)

## checking correlations
round(cor(dat[,-c(1)], use = "pairwise"), 3)

## remove colliniar variables
dat$gai <- NULL
dat$gbgr <- NULL
dat$gai_r <- NULL
dat$gbgr_r <- NULL
dat$perc_stigma <- NULL
dat$age <- NULL
dat$grpid <- NULL
dat$perc_stigma_gm <- NULL
dat$age_gm <- NULL
dat$grpid_gm <- NULL
dat$polori <- NULL
dat$polori_gm <- NULL
dat$gini2016 <- NULL
dat$gdp2016 <- NULL


#############################
## multiple imputation step
#############################

dat$gen_min <- as.factor(dat$gen_min)

set.seed(111) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen_min'), idvars = c('Sample'), m = 5)
summary(impute.out)
 
a.mids <- datlist2mids(impute.out$imputations)

str(impute.out)

########################################
## Main Bayesian regression analysis
########################################

####
## Individual Level
####

## Setting proper priors fro misspecified priors in the default configuration.

mypriors_perc <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 10), class = "b", coef = "gen_min1"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              #prior(normal(0, 10), class = "b", coef = "polori_gmc"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )

mypriors_ins <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 10), class = "b", coef = "gen_min1"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              #prior(normal(0, 10), class = "b", coef = "polori_gmc"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )

wilac_perc <- brm_multiple(wilac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                               ## control variables
                               gen_min +
                               age_gmc +
                               grpid_gmc +
                               #polori_gmc +
                               gini2016s +
                               gdp2016s +
                               (1 | Sample),
                           data = a.mids,
                           prior = mypriors_perc,
                           cores = 3,
                           seed = 123,
                           refresh = 0,
                           open_progress = FALSE,
                           save_pars = save_pars(all = TRUE),
                                        #file = 'models/wilac_perc'
                           )
summary(wilac_perc)

##Participation by perceive stigma model
colac_perc <- brm_multiple(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              gen_min +
                              age_gmc +
                              grpid_gmc +
                              #polori_gmc +
                              gini2016s +
                              gdp2016s +
                              (1 | Sample),
                          data = a.mids,
                          prior = mypriors_perc,
                          cores = 3,
                          seed = 222,
                          refresh = 0,
                          open_progress = FALSE,
                          save_pars = save_pars(all = TRUE),
                          #file = 'models/colac_perc'
                          )
summary(colac_perc)

##Willingness by institutional stigma model
wilac_ins <- brm_multiple(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             gen_min +
                             age_gmc +
                             grpid_gmc +
                             #polori_gmc +
                             gini2016s +
                             gdp2016s +
                             (1 | Sample),
                         data = a.mids,
                         prior = mypriors_ins,
                         cores = 3,
                         seed = 111,
                         refresh = 0,
                         open_progress = FALSE,
                         save_pars = save_pars(all = TRUE),
                         # file = 'models/wilac_ins'
                         )
summary(wilac_ins)

##Participation by institutional stigma model

colac_ins <- brm_multiple(colac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             gen_min +
                             age_gmc +
                             grpid_gmc +
                             #polori_gmc +
                             gini2016s +
                             gdp2016s +
                             (1 | Sample),
                         data = a.mids,
                         prior = mypriors_ins,
                         cores = 3,
                         seed = 111,
                         refresh = 0,
                         open_progress = FALSE,
                         save_pars = save_pars(all = TRUE),
                         # file = 'models/colac_ins'
                         )
summary(colac_ins)

########################################
## Main Bayesian regression analysis with Political Orientation
########################################

####
## Individual Level
####

## Setting proper priors fro misspecified priors in the default configuration.

mypriors_perc <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 10), class = "b", coef = "gen_min1"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              prior(normal(0, 10), class = "b", coef = "polori_gmc"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )

mypriors_ins <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 10), class = "b", coef = "gen_min1"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              prior(normal(0, 10), class = "b", coef = "polori_gmc"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )

wilac_perc_po <- brm_multiple(wilac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                               ## control variables
                               gen_min +
                               age_gmc +
                               grpid_gmc +
                               polori_gmc +
                               gini2016s +
                               gdp2016s +
                               (1 | Sample),
                           data = a.mids,
                           prior = mypriors_perc,
                           cores = 3,
                           seed = 123,
                           refresh = 0,
                           open_progress = FALSE,
                           save_pars = save_pars(all = TRUE),
                                        #file = 'models/wilac_perc'
                           )
summary(wilac_perc_po)

##Participation by perceive stigma model
colac_perc_po <- brm_multiple(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              gen_min +
                              age_gmc +
                              grpid_gmc +
                              polori_gmc +
                              gini2016s +
                              gdp2016s +
                              (1 | Sample),
                          data = a.mids,
                          prior = mypriors_perc,
                          cores = 3,
                          seed = 222,
                          refresh = 0,
                          open_progress = FALSE,
                          save_pars = save_pars(all = TRUE),
                          #file = 'models/colac_perc'
                          )
summary(colac_perc_po)

##Willingness by institutional stigma model
wilac_ins_po <- brm_multiple(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             gen_min +
                             age_gmc +
                             grpid_gmc +
                             polori_gmc +
                             gini2016s +
                             gdp2016s +
                             (1 | Sample),
                         data = a.mids,
                         prior = mypriors_ins,
                         cores = 3,
                         seed = 111,
                         refresh = 0,
                         open_progress = FALSE,
                         save_pars = save_pars(all = TRUE),
                         # file = 'models/wilac_ins'
                         )
summary(wilac_ins_po)

##Participation by institutional stigma model

colac_ins_po <- brm_multiple(colac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             gen_min +
                             age_gmc +
                             grpid_gmc +
                             polori_gmc +
                             gini2016s +
                             gdp2016s +
                             (1 | Sample),
                         data = a.mids,
                         prior = mypriors_ins,
                         cores = 3,
                         seed = 111,
                         refresh = 0,
                         open_progress = FALSE,
                         save_pars = save_pars(all = TRUE),
                         # file = 'models/colac_ins'
                         )
summary(colac_ins_po)

wilac_perc <- add_criterion(wilac_perc, "loo")
wilac_perc_po <- add_criterion(wilac_perc_po, "loo")
loo_compare(wilac_perc, wilac_perc_po, criterion = "loo")

wilac_perc <- add_criterion(wilac_perc, "waic")
wilac_perc_po <- add_criterion(wilac_perc_po, "waic")
loo_compare(wilac_perc, wilac_perc_po, criterion = "waic")
