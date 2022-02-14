library(ggplot2)
library(Amelia)
library(sjPlot)
library(psych)
library(gridExtra)
library(mice)
library(miceadds)
library(brms)
library(ggrepel)

dat <- read.csv('data/fulldataset.csv', row.names = 1)
names(dat)
str(dat)

dat$polori[dat$polori == 77] <- NA
dat$polori[dat$polori == 66] <- NA
dat$polori[dat$polori == 8] <- NA

table(dat$polori)

prop.table(table(dat$sexo))

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

plot_between <- function(model, x, y, xlab, ylab, data_within, data_between){
    g <- plot(conditional_effects(model, x), plot = FALSE)[[1]] +
        geom_point(
            aes_string(x = x, y = y), 
            data = data_within,
            size = 4,
            alpha = .05,
            color = dot_colour,
            stroke = 0,
            position = position_jitter(w = 0.15, h = 0.15),
            inherit.aes = FALSE) +
        geom_point(
            aes_string(x = x, y = y), 
            data = data_between,
            size = 4,
            color = country_dot_colour,
            stroke = 0,
            inherit.aes = FALSE) +
        geom_text_repel(
            aes_string(x = x, y = y, label = "Sample"),
            max.overlaps = Inf,
            data = data_between,
            inherit.aes = FALSE) +
        theme_bw() +
        ylim(1, 7) +
        xlab(xlab) +
        ylab(ylab)
    return(g)
}

plot_within <- function(model, x, y, xlab, ylab, data_within){
    g <- plot(conditional_effects(model, x), plot = FALSE)[[1]] +
        geom_point(
            aes_string(x = x, y = y), 
            data = data_within,
            size = 4,
            alpha = .05,
            color = dot_colour,
            stroke = 0,
            position = position_jitter(w = 0.15, h = 0.15),
            inherit.aes = FALSE) +
        theme_bw() +
        ylim(1, 7) +
        xlab(xlab) +
        ylab(ylab)
    return(g)
}

#################
## descriptives
#################

options(browser="/usr/bin/firefox")

tab_df(aggregate(ID ~ Sample, data = dat, FUN = length))
prop.table(table(dat$gen_min))*100
prop.table(table(dat$new_gen))*100

frqdat <- subset(dat, select = -c(ID, Sample, gen_min))
tab_df(describe(frqdat), digits = 2)
tab_corr(frqdat,  digits = 2)

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

#############################
## multiple imputation step
#############################

dat$gen_min <- as.factor(dat$gen_min)
dat$new_gen <- as.factor(dat$new_gen)
set.seed(111) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen_min', 'new_gen'), idvars = c('Sample'), m = 5)
summary(impute.out)
 
a.mids <- datlist2mids(impute.out$imputations)

str(impute.out)

########################################
## Main Bayesian regression analysis
########################################

#################
## Within models
#################

## Setting proper priors fro misspecified priors in the default configuration.
## Political orientation presented problems that prevented convergence of the models.\
## Those problems where unsolvable and thus we dropped the variable.

mypriors_perc <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 10), class = "b", coef = "gen_min1"),
              prior(normal(0, 10), class = "b", coef = "new_gen2"),
              prior(normal(0, 10), class = "b", coef = "new_gen3"),
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
              prior(normal(0, 10), class = "b", coef = "new_gen2"),
              prior(normal(0, 10), class = "b", coef = "new_gen3"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              #prior(normal(0, 10), class = "b", coef = "polori_gmc"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )

##Intention of Participation by perceive stigma model
wilac_perc <- brm_multiple(wilac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                               ## control variables
                               gen_min +
                               new_gen +
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
                           file = 'models/wilac_perc'
                           )

summary(wilac_perc)
bayestestR::hdi(wilac_perc, ci = c(.90, .95))
bayes_R2(wilac_perc)

##Participation by perceive stigma model
colac_perc <- brm_multiple(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              gen_min +
                              new_gen +
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
                          file = 'models/colac_perc'
                          )

summary(colac_perc)
bayestestR::hdi(colac_perc, ci = c(.90, .95))
bayes_R2(colac_perc)

##Linear model for comparison

mypriors_perc_linear <- c(prior(normal(0, 10), class = "Intercept"),
                          prior(normal(0, 10), class = "b", coef = "perc_stigma_gmc"),
                          prior(normal(0, 10), class = "b", coef = "age_gmc"),
                          prior(normal(0, 10), class = "b", coef = "gen_min1"),
                          prior(normal(0, 10), class = "b", coef = "new_gen2"),
                          prior(normal(0, 10), class = "b", coef = "new_gen3"),
                          prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
                          #prior(normal(0, 10), class = "b", coef = "polori_gmc"),
                          prior(normal(0, 10), class = "b", coef = "gdp2016s"),
                          prior(normal(0, 10), class = "b", coef = "gini2016s")
                          )

mypriors_ins_linear <- c(prior(normal(0, 10), class = "Intercept"),
                         prior(normal(0, 10), class = "b", coef = "ins_stigma"),
                         prior(normal(0, 10), class = "b", coef = "age_gmc"),
                         prior(normal(0, 10), class = "b", coef = "gen_min1"),
                         prior(normal(0, 10), class = "b", coef = "new_gen2"),
                         prior(normal(0, 10), class = "b", coef = "new_gen3"),
                         prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
                         #prior(normal(0, 10), class = "b", coef = "polori_gmc"),
                         prior(normal(0, 10), class = "b", coef = "gdp2016s"),
                         prior(normal(0, 10), class = "b", coef = "gini2016s")
                         )

wilac_perc_linear <- brm_multiple(wilac ~ perc_stigma_gmc +
                                      ## control variables
                                      gen_min +
                                      new_gen +
                                      age_gmc +
                                      grpid_gmc +
                                      #polori_gmc +
                                      gini2016s +
                                      gdp2016s +
                                      (1 | Sample),
                                  data = a.mids,
                                  prior = mypriors_perc_linear,
                                  cores = 3,
                                  seed = 111,
                                  refresh = 0,
                                  open_progress = FALSE,
                                  save_pars = save_pars(all = TRUE),
                                  file = 'models/wilac_perc_linear'
                                  )

colac_perc_linear <- brm_multiple(colac ~ perc_stigma_gmc +
                                      ## control variables
                                      gen_min +
                                      new_gen +
                                      age_gmc +
                                      grpid_gmc +
                                      #polori_gmc +
                                      gini2016s +
                                      gdp2016s +
                                      (1 | Sample),
                                  data = a.mids,
                                  prior = mypriors_perc_linear,
                                  cores = 3,
                                  seed = 111,
                                  refresh = 0,
                                  open_progress = FALSE,
                                  save_pars = save_pars(all = TRUE),
                                  file = 'models/colac_perc_linear'
                                  )


bayes_factor(wilac_perc, wilac_perc_linear)
bayes_factor(colac_perc, colac_perc_linear)

##################
## Between Models
##################

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

mypriors_ins_g <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE2"),
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
summary(wilac_perc_g)
bayestestR::hdi(wilac_perc_g, ci = c(.90, .95))
bayes_R2(wilac_perc_g)

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
summary(colac_perc_g)
bayestestR::hdi(colac_perc_g, ci = c(.90, .95))
bayes_R2(colac_perc_g)

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
                            prior = mypriors_ins_g,
                            cores = 3,
                            seed = 123,
                            refresh = 0,
                            open_progress = FALSE,
                            save_pars = save_pars(all = TRUE),
                            file = 'models/wilac_ins_g'
                            )
summary(wilac_ins_g)
bayestestR::hdi(wilac_ins_g, ci = c(.89, .90, .95))
bayes_R2(wilac_ins_g)

colac_ins_g <- brm_multiple(colac ~ poly(ins_stigma, 2, raw=FALSE) +
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
                            prior = mypriors_ins_g,
                            cores = 3,
                            seed = 123,
                            refresh = 0,
                            open_progress = FALSE,
                            save_pars = save_pars(all = TRUE),
                            file = 'models/colac_ins_g'
                            )
summary(colac_ins_g)
bayestestR::hdi(colac_ins_g, ci = c(.90, .95))
bayes_R2(colac_ins_g)


## Linear models for comparison

mypriors_perc_g_l <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "perc_stigma_gm"),
              prior(normal(0, 10), class = "b", coef = "perc_stigma_gmc"),
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

mypriors_ins_g_l <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "ins_stigma"),
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


wilac_perc_g_l <- brm_multiple(wilac ~ perc_stigma_gm +
                                 perc_stigma_gmc +
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
                             prior = mypriors_perc_g_l,
                             cores = 3,
                             seed = 123,
                             refresh = 0,
                             open_progress = FALSE,
                             save_pars = save_pars(all = TRUE),
                             file = 'models/wilac_perc_g_l'
                             )
summary(wilac_perc_g_l)

colac_perc_g_l <- brm_multiple(colac ~ perc_stigma_gm +
                                 perc_stigma_gmc +
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
                             prior = mypriors_perc_g_l,
                             cores = 3,
                             seed = 123,
                             refresh = 0,
                             open_progress = FALSE,
                             save_pars = save_pars(all = TRUE),
                             file = 'models/colac_perc_g_l'
                             )
summary(colac_perc_g_l)

wilac_ins_g_l <- brm_multiple(wilac ~ ins_stigma +
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
                            prior = mypriors_ins_g_l,
                            cores = 3,
                            seed = 123,
                            refresh = 0,
                            open_progress = FALSE,
                            save_pars = save_pars(all = TRUE),
                            file = 'models/wilac_ins_g_l'
                            )
summary(wilac_ins_g_l)


colac_ins_g_l <- brm_multiple(colac ~ ins_stigma +
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
                            prior = mypriors_ins_g_l,
                            cores = 3,
                            seed = 123,
                            refresh = 0,
                            open_progress = FALSE,
                            save_pars = save_pars(all = TRUE),
                            file = 'models/colac_ins_g_l'
                            )
summary(colac_ins_g_l)

bayes_factor(wilac_perc_g, wilac_perc_g_l)
bayes_factor(colac_perc_g, colac_perc_g_l)
bayes_factor(wilac_ins_g, wilac_ins_g_l)
bayes_factor(colac_ins_g, colac_ins_g_l)


##########
## Plots
##########

dot_colour <- "black"
country_dot_colour <- "tomato3"

## Within Plots

g1 <- plot_within(model = wilac_perc,
                  x = "perc_stigma_gmc",
                  y = "wilac",
                  xlab = "Perceived Stigma",
                  ylab = "Collective Action Intentions",
                  data_within = dat)

g2 <- plot_within(model = colac_perc,
                  x = "perc_stigma_gmc",
                  y = "colac",
                  xlab = "Perceived Stigma",
                  ylab = "Participation in Collective Action",
                  data_within = dat)

ggsave("plots/perc_within.png", arrangeGrob(g1, g2, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)

## Between Plots

country_means <- aggregate(cbind(perc_stigma_gm, ins_stigma, wilac, colac) ~ Sample,
                           data = dat,
                           FUN = mean)

g3 <- plot_between(model = wilac_perc_g,
                   x = "perc_stigma_gm",
                   y = "wilac",
                   xlab = "Perceived Stigma",
                   ylab = "Collective Action Intentions",
                   data_within = dat,
                   data_between = country_means)
g4 <- plot_between(model = colac_perc_g,
                   x = "perc_stigma_gm",
                   y = "colac",
                   xlab = "Perceived Stigma",
                   ylab = "Participation in Collective Action",
                   data_within = dat,
                   data_between = country_means)
g5 <- plot_between(model = wilac_ins_g,
                   x = "ins_stigma",
                   y = "wilac",
                   xlab = "Institutional Stigma",
                   ylab = "Collective Action Intentions",
                   data_within = dat,
                   data_between = country_means)
g6 <- plot_between(model = colac_ins_g,
                   x = "ins_stigma",
                   y = "colac",
                   xlab = "Institutional Stigma",
                   ylab = "Participation in Collective Action",
                   data_within = dat,
                   data_between = country_means)

ggsave("plots/perc_between.png", arrangeGrob(g3, g4, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)
ggsave("plots/ins_between.png", arrangeGrob(g5, g6, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)


######################
##Tipping point values
######################
get_tipping_point(g3)
get_tipping_point(g4)
get_tipping_point(g5)
get_tipping_point(g6)

#######################
## Mediation analysis
#######################

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

## remove colliniar variables
dat$gai <- NULL
dat$gbgr <- NULL
dat$gai_r <- NULL
dat$gbgr_r <- NULL
dat$perc_stigma <- NULL
dat$perc_stigma_gmc <- NULL
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

#############################
## multiple imputation step
#############################

dat$gen_min <- as.factor(dat$gen_min)
dat$new_gen <- as.factor(dat$new_gen)
set.seed(111) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen_min', 'new_gen'), idvars = c('Sample'), m = 5)
summary(impute.out)
 
a.mids <- datlist2mids(impute.out$imputations)


mx_priors <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "ins_stigma"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 10), class = "b", coef = "gen_min1"),
              prior(normal(0, 10), class = "b", coef = "new_gen2"),
              prior(normal(0, 10), class = "b", coef = "new_gen3"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              #prior(normal(0, 10), class = "b", coef = "polori_gmc"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )

ymx_priors <- c(prior(normal(0, 10), class = "Intercept"),
                prior(normal(0, 10), class = "b", coef = "polyperc_stigma_mc2rawEQFALSE1"),
                prior(normal(0, 10), class = "b", coef = "polyperc_stigma_mc2rawEQFALSE2"),
                prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE1"),
                prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE2"),
                prior(normal(0, 10), class = "b", coef = "age_gmc"),
                prior(normal(0, 10), class = "b", coef = "gen_min1"),
                prior(normal(0, 10), class = "b", coef = "new_gen2"),
                prior(normal(0, 10), class = "b", coef = "new_gen3"),
                prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
                #prior(normal(0, 10), class = "b", coef = "polori_gmc"),
                prior(normal(0, 10), class = "b", coef = "gdp2016s"),
                prior(normal(0, 10), class = "b", coef = "gini2016s")
                )

mxmodel <- brm_multiple(perc_stigma_mc ~ ins_stigma +
                            ## control variables
                            gen_min +
                            new_gen +
                            age_gmc +
                            grpid_gmc +
                            #polori_gmc +
                            gini2016s +
                            gdp2016s +
                            (1 | Sample),
                        data = a.mids,
                        prior = mx_priors,
                        cores = 2,
                        file = 'models/mxmodel')
summary(mxmodel)
bayestestR::hdi(mxmodel, ci = c(.90, .95))
bayes_R2(mxmodel)

ymxmodel_wi <-  brm_multiple(wilac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                                 poly(ins_stigma, 2, raw = FALSE)  +
                                 ## control variables
                                 gen_min +
                                 new_gen +
                                 age_gmc +
                                 grpid_gmc +
                                 #polori_gmc +
                                 gini2016s +
                                 gdp2016s +
                                 (1 | Sample),
                             data = a.mids,
                             prior = ymx_priors,
                             cores = 2,
                             file = 'models/ymxmodel_wi')
summary(ymxmodel_wi)
bayestestR::hdi(ymxmodel_wi, ci = c(.90, .95))
bayes_R2(ymxmodel_wi)

ymxmodel_co <- brm_multiple(colac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                                poly(ins_stigma, 2, raw = FALSE)  +
                                ## control variables
                                gen_min +
                                new_gen +
                                age_gmc +
                                grpid_gmc +
                                #polori_gmc +
                                gini2016s +
                                gdp2016s +
                                (1 | Sample),
                            data = a.mids,
                            prior = ymx_priors,
                            cores = 2,
                            file = 'models/ymxmodel_co')
summary(ymxmodel_co)
bayestestR::hdi(ymxmodel_co, ci = c(.90, .95))
bayes_R2(ymxmodel_co)

###############################################
## instantaneous indirect effects calculation
###############################################

####
## Participation in collective actions
####

a <- as.vector(posterior_samples(mxmodel, 'ins_stigma', as.matrix = T))
i1 <- as.vector(posterior_samples(mxmodel, 'b_Intercept', as.matrix = T))
b1 <- as.vector(posterior_samples(ymxmodel_co, 'polyperc_stigma_mc2rawEQFALSE1', as.matrix=T))
b2 <- as.vector(posterior_samples(ymxmodel_co, 'polyperc_stigma_mc2rawEQFALSE2', as.matrix=T))

df <- data.frame()
means_co <- c()
xvals <- seq(-1.17764, 3.22305, length.out = 100)
for (val in xvals){
    b <- b1 + 2*b2*(i1 + a*val + sum(mean(dat$age_gmc, na.rm=T),
                                     mean(dat$grpid_gmc, na.rm=T),
                                     mean(dat$polori_gmc, na.rm=T),
                                     mean(dat$gini2016s, na.rm=T),
                                     mean(dat$gdp2016s, na.rm=T)))
    theta <- a*b
    out_ith <- quantile(theta, probs = c(0.025, 0.05, 0.5, 0.95, 0.975), names=F)
    means_co <- c(means_co, mean(out_ith))
    df <- rbind(df, out_ith)
}
names <- c('CIll', 'CIl', 'Median', 'CIu', 'CIuu')
colnames(df) <- names

medplot_co <- ggplot(data = df) +
    geom_line(aes(x=xvals, y=means_co)) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
    geom_ribbon(aes(x=xvals, ymin=CIl, ymax=CIu),
                alpha=0.3,
                fill='blue') +
    geom_ribbon(aes(x=xvals, ymin=CIll, ymax=CIuu),
                alpha=0.3,
                fill='red') +
    xlab('Institutional Stigma') +
    ylab('Instantaneous Indirect Effect on Participation in Collective Action') +
    theme_bw() +
    ylim(-20, 20)

####
## Collective actions intentions
####

a <- as.vector(posterior_samples(mxmodel, 'ins_stigma', as.matrix = T))
i1 <- as.vector(posterior_samples(mxmodel, 'b_Intercept', as.matrix = T))
b1 <- as.vector(posterior_samples(ymxmodel_wi, 'polyperc_stigma_mc2rawEQFALSE1', as.matrix=T))
b2 <- as.vector(posterior_samples(ymxmodel_wi, 'polyperc_stigma_mc2rawEQFALSE2', as.matrix=T))

df <- data.frame()
means_wi <- c()
xvals <- seq(-1.17764, 3.22305, length.out = 100)
for (val in xvals){
    b <- b1 + 2*b2*(i1 + a*val + sum(mean(dat$age_gmc, na.rm=T),
                                     mean(dat$grpid_gmc, na.rm=T),
                                     mean(dat$polori_gmc, na.rm=T),
                                     mean(dat$gini2016s, na.rm=T),
                                     mean(dat$gdp2016s, na.rm=T)))
    theta <- a*b
    out_ith <- quantile(theta, probs = c(0.025, 0.05, 0.5, 0.95, 0.975), names=F)
    means_wi <- c(means_wi, mean(out_ith))
    df <- rbind(df, out_ith)
}
names <- c('CIll', 'CIl', 'Median', 'CIu', 'CIuu')
colnames(df) <- names

medplot_wi <- ggplot(data = df) +
    geom_line(aes(x=xvals, y=means_wi)) +
    geom_abline(slope = 0, intercept = 0, linetype = "dashed") +
    geom_ribbon(aes(x=xvals, ymin=CIl, ymax=CIu),
                alpha=0.3,
                fill='blue') +
    geom_ribbon(aes(x=xvals, ymin=CIll, ymax=CIuu),
                alpha=0.3,
                fill='red') +
    xlab('Institutional Stigma') +
    ylab('Instantaneous Indirect Effect on Collective Action Intentions') +
    theme_bw() +
    ylim(-20, 20)

######################

ggsave("plots/medplot.png", arrangeGrob(medplot_wi, medplot_co, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)
