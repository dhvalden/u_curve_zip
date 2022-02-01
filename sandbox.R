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
             'new_gen',
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

## auxiliary functions

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

get_tipping_point <- function(plot_object){
    plot_data = ggplot_build(plot_object)$data[[1]]
    x_max_y_value = plot_data$x[which.max(plot_data$y)]
    return(x_max_y_value)
}

## creating centered variables
dat <- get_mean_centered('perc_stigma', 'Sample', dat)
dat <- get_mean_centered('age', 'Sample', dat)
dat <- get_mean_centered('polori', 'Sample', dat)
dat <- get_mean_centered('grpid', 'Sample', dat)
dat$ins_stigma <- as.vector(scale(dat$ins_stigma))
dat$gini2016s <- as.vector(scale(dat$gini2016))
dat$gdp2016s <- as.vector(scale(dat$gdp2016))

Gen_minProp <- prop.table(table(dat$Sample, dat$gen_min), 1)
Gen_minProp <- as.data.frame.matrix(Gen_minProp)
colnames(Gen_minProp) <- c('prop_no_gen_min', 'prop_gen_min')
Gen_minProp$Sample <- row.names(Gen_minProp)

new_gen_prop <- prop.table(table(dat$Sample, dat$new_gen), 1)
new_gen_prop <- as.data.frame.matrix(new_gen_prop)
colnames(new_gen_prop) <- c('prop_male', 'prop_female', 'prop_other')
new_gen_prop$Sample <- row.names(new_gen_prop)

dat <- merge(dat, Gen_minProp, by='Sample')
dat <- merge(dat, new_gen_prop, by='Sample')


str(dat)
names(dat)
summary(dat)

## checking correlations
round(cor(dat[,-c(1, 2, 3, 4)], use = "pairwise"), 3)

## remove colliniar variables
dat$gai <- NULL
dat$gbgr <- NULL
dat$gai_r <- NULL
dat$gbgr_r <- NULL
dat$perc_stigma <- NULL
dat$perc_stigma_mc <- NULL
dat$age <- NULL
dat$grpid <- NULL
dat$age_mc <- NULL
dat$grpid_mc <- NULL
dat$polori <- NULL
dat$polori_mc <- NULL
dat$gini2016 <- NULL
dat$gdp2016 <- NULL
dat$prop_no_gen_min <- NULL
dat$prop_male <- NULL


str(dat)

round(cor(dat[,-c(1, 2)], use = "pairwise"), 3)

#############################
## multiple imputation step
#############################

dat$gen_min <- as.factor(dat$gen_min)
dat$new_gen <- as.factor(dat$new_gen)
set.seed(111) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen_min', 'new_gen'), idvars = c('Sample'), m = 5)
summary(impute.out)
 
a.mids <- datlist2mids(impute.out$imputations)


mypriors_perc_g <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gm2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gm2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "gen_min1"),
              prior(normal(0, 10), class = "b", coef = "prop_gen_min"),
              prior(normal(0, 10), class = "b", coef = "new_gen2"),
              prior(normal(0, 10), class = "b", coef = "new_gen3"),
              prior(normal(0, 10), class = "b", coef = "prop_female"),
              prior(normal(0, 10), class = "b", coef = "prop_other"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 10), class = "b", coef = "age_gm"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              prior(normal(0, 10), class = "b", coef = "grpid_gm"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )


wilac_perc_g <- brm_multiple(wilac ~ poly(perc_stigma_gm, 2, raw=FALSE) +
                                 poly(perc_stigma_gmc, 2, raw=FALSE) +
                                 gen_min +
                                 prop_gen_min +
                                 new_gen +
                                 prop_female +
                                 prop_other +
                                 age_gmc +
                                 age_gm +
                                 grpid_gmc +
                                 grpid_gm +
                                 gini2016s +
                                 gdp2016s +
                                 (1 | Sample),
                             data = a.mids,
                             prior = mypriors_perc_g,
                             cores = 3,
                             seed = 123,
                             refresh = 0,
                             open_progress = FALSE,
                             save_pars = save_pars(all = TRUE),
                             file = 'models/wilac_perc_g'
                             )

summary(wilac_perc_g, prob = .9)

colac_perc_g <- brm_multiple(colac ~ poly(perc_stigma_gm, 2, raw=FALSE) +
                                 poly(perc_stigma_gmc, 2, raw=FALSE) +
                                 gen_min +
                                 prop_gen_min +
                                 new_gen +
                                 prop_female +
                                 prop_other +
                                 age_gmc +
                                 age_gm +
                                 grpid_gmc +
                                 grpid_gm +
                                 gini2016s +
                                 gdp2016s +
                                 (1 | Sample),
                             data = a.mids,
                             prior = mypriors_perc_g,
                             cores = 3,
                             seed = 123,
                             refresh = 0,
                             open_progress = FALSE,
                             save_pars = save_pars(all = TRUE),
                             file = 'models/colac_perc_g'
                             )
summary(colac_perc_g, prob = .9)

wilac_ins_g <- brm_multiple(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                                 gen_min +
                                 prop_gen_min +
                                 new_gen +
                                 prop_female +
                                 prop_other +
                                 age_gmc +
                                 age_gm +
                                 grpid_gmc +
                                 grpid_gm +
                                 gini2016s +
                                 gdp2016s +
                                 (1 | Sample),
                             data = a.mids,
                             cores = 3,
                             seed = 123,
                             refresh = 0,
                             open_progress = FALSE,
                             save_pars = save_pars(all = TRUE),
                             )
summary(wilac_ins_g, prob = .9)
