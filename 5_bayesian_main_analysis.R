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
             'gen',
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

options(browser="/usr/bin/google-chrome-stable")

tab_df(aggregate(ID ~ Sample, data = dat, FUN = length))
prop.table(table(dat$gen))*100

frqdat <- subset(dat, select = -c(ID, Sample, gen))
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

dat$gen <- as.factor(dat$gen)

set.seed(111) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen'), idvars = c('Sample'), m = 5)
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
              prior(normal(0, 20), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE1"),
              prior(normal(0, 20), class = "b", coef = "polyperc_stigma_gmc2rawEQFALSE2"),
              prior(normal(0, 20), class = "b", coef = "age_gmc"),
              prior(normal(0, 2), class = "b", coef = "gen2"),
              prior(normal(0, 2), class = "b", coef = "gen3"),
              prior(normal(0, 2), class = "b", coef = "gen4"),
              prior(normal(0, 2), class = "b", coef = "gen5"),
              prior(normal(0, 20), class = "b", coef = "grpid_gmc"),
              prior(normal(0, 20), class = "b", coef = "polori_gmc"),
              prior(normal(0, 20), class = "b", coef = "gdp2016s"),
              prior(normal(0, 20), class = "b", coef = "gini2016s")
              )

mypriors_ins <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 20), class = "b", coef = "polyins_stigma2rawEQFALSE1"),
              prior(normal(0, 20), class = "b", coef = "polyins_stigma2rawEQFALSE2"),
              prior(normal(0, 20), class = "b", coef = "age_gmc"),
              prior(normal(0, 2), class = "b", coef = "gen2"),
              prior(normal(0, 2), class = "b", coef = "gen3"),
              prior(normal(0, 2), class = "b", coef = "gen4"),
              prior(normal(0, 2), class = "b", coef = "gen5"),
              prior(normal(0, 20), class = "b", coef = "grpid_gmc"),
              prior(normal(0, 20), class = "b", coef = "polori_gmc"),
              prior(normal(0, 20), class = "b", coef = "gdp2016s"),
              prior(normal(0, 20), class = "b", coef = "gini2016s")
              )


wilac_perc <- brm_multiple(wilac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              gen +
                              age_gmc +
                              grpid_gmc +
                              polori_gmc +
                              gini2016s +
                              gdp2016s +
                              (1 | Sample),
                          data = a.mids,
                          prior = mypriors_perc,
                          cores = 2,
                          save_pars = save_pars(all = TRUE),
                          silent = 2,
                          file = 'models/wilac_perc')
summary(wilac_perc)

##Participation by perceive stigma model
colac_perc <- brm_multiple(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              gen +
                              age_gmc +
                              grpid_gmc +
                              polori_gmc +
                              gini2016s +
                              gdp2016s +
                              (1 | Sample),
                          data = a.mids,
                          prior = mypriors_perc,
                          cores = 2,
                          save_pars = save_pars(all = TRUE),
                          silent = 2,
                          file = 'models/colac_perc')
summary(colac_perc)

##Willingness by institutional stigma model
wilac_ins <- brm_multiple(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             gen +
                             age_gmc +
                             grpid_gmc +
                             polori_gmc +
                             gini2016s +
                             gdp2016s +
                             (1 | Sample),
                         data = a.mids,
                         prior = mypriors_ins,
                         cores = 2,
                         save_pars = save_pars(all = TRUE),
                         silent = 2,
                         file = 'models/wilac_ins')
summary(wilac_ins)

##Participation by institutional stigma model

colac_ins <- brm_multiple(colac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             gen +
                             age_gmc +
                             grpid_gmc +
                             polori_gmc +
                             gini2016s +
                             gdp2016s +
                             (1 | Sample),
                         data = a.mids,
                         prior = mypriors_ins,
                         cores = 2,
                         save_pars = save_pars(all = TRUE),
                         silent = 2,
                         file = 'models/colac_ins')
summary(colac_ins)

tab_model(wilac_perc, colac_perc, show.se = TRUE)
tab_model(wilac_ins, colac_ins, show.se = TRUE)

##Linear model for comparison

mypriors_perc_linear <- c(prior(normal(0, 10), class = "Intercept"),
                          prior(normal(0, 20), class = "b", coef = "perc_stigma_gmc"),
                          prior(normal(0, 20), class = "b", coef = "age_gmc"),
                          prior(normal(0, 2), class = "b", coef = "gen2"),
                          prior(normal(0, 2), class = "b", coef = "gen3"),
                          prior(normal(0, 2), class = "b", coef = "gen4"),
                          prior(normal(0, 2), class = "b", coef = "gen5"),
                          prior(normal(0, 20), class = "b", coef = "grpid_gmc"),
                          prior(normal(0, 20), class = "b", coef = "polori_gmc"),
                          prior(normal(0, 20), class = "b", coef = "gdp2016s"),
                          prior(normal(0, 20), class = "b", coef = "gini2016s")
                          )

mypriors_ins_linear <- c(prior(normal(0, 10), class = "Intercept"),
                         prior(normal(0, 20), class = "b", coef = "ins_stigma"),
                         prior(normal(0, 20), class = "b", coef = "age_gmc"),
                         prior(normal(0, 2), class = "b", coef = "gen2"),
                         prior(normal(0, 2), class = "b", coef = "gen3"),
                         prior(normal(0, 2), class = "b", coef = "gen4"),
                         prior(normal(0, 2), class = "b", coef = "gen5"),
                         prior(normal(0, 20), class = "b", coef = "grpid_gmc"),
                         prior(normal(0, 20), class = "b", coef = "polori_gmc"),
                         prior(normal(0, 20), class = "b", coef = "gdp2016s"),
                         prior(normal(0, 20), class = "b", coef = "gini2016s")
                         )

wilac_perc_linear <- brm_multiple(wilac ~ perc_stigma_gmc +
                                      ## control variables
                                      gen +
                                      age_gmc +
                                      grpid_gmc +
                                      polori_gmc +
                                      gini2016s +
                                      gdp2016s +
                                      (1 | Sample),
                                  data = a.mids,
                                  prior = mypriors_perc_linear,
                                  cores = 2,
                                  save_pars = save_pars(all = TRUE),
                                  silent = 2,
                                  file = 'models/wilac_perc_linear')

colac_perc_linear <- brm_multiple(colac ~ perc_stigma_gmc +
                                      ## control variables
                                      gen +
                                      age_gmc +
                                      grpid_gmc +
                                      polori_gmc +
                                      gini2016s +
                                      gdp2016s +
                                      (1 | Sample),
                                  data = a.mids,
                                  prior = mypriors_perc_linear,
                                  cores = 2,
                                  save_pars = save_pars(all = TRUE),
                                  silent = 2,
                                  file = 'models/colac_perc_linear')

wilac_ins_linear <- brm_multiple(wilac ~ ins_stigma +
                                     ## control variables
                                     gen +
                                     age_gmc +
                                     grpid_gmc +
                                     polori_gmc +
                                     gini2016s +
                                     gdp2016s +
                                     (1 | Sample),
                                 data = a.mids,
                                 prior = mypriors_ins_linear,
                                 cores = 2,
                                 save_pars = save_pars(all = TRUE),
                                 silent = 2,
                                 file = 'models/wilac_ins_linear')

colac_ins_linear <- brm_multiple(colac ~ ins_stigma +
                                     ## control variables
                                     gen +
                                     age_gmc +
                                     grpid_gmc +
                                     polori_gmc +
                                     gini2016s +
                                     gdp2016s +
                                     (1 | Sample),
                                 data = a.mids,
                                 prior = mypriors_ins_linear,
                                 cores = 2,
                                 save_pars = save_pars(all = TRUE),
                                 silent = 2,
                                 file = 'models/colac_ins_linear')


bayes_factor(wilac_perc, wilac_perc_linear)
bayes_factor(colac_perc, colac_perc_linear)

bayes_factor(wilac_ins, wilac_ins_linear)
bayes_factor(colac_ins, colac_ins_linear)


####
##Country Level
###

## getting complete data

completedData <- complete(a.mids, 'long')
dat_imp <- aggregate(completedData[-c(1,2)], by=list(completedData$ID), FUN=mean)
dat_imp$Group.1 <- NULL
dat_imp$Sample <- NULL
dat_imp$gen <- NULL

dat_imp <- merge(dat[c('ID','Sample','gen')], dat_imp, by='ID')

wilac <- aggregate(wilac ~ Sample, data = dat_imp, mean)
colnames(wilac) <- c('Sample', 'MWilac')
colac <- aggregate(colac ~ Sample, data = dat_imp, mean)
colnames(colac) <- c('Sample', 'MColac')
age <- aggregate(age_mc ~ Sample, data = dat_imp, mean)
colnames(age) <- c('Sample', 'MAge')
grpid <- aggregate(grpid_mc ~ Sample, data = dat_imp, mean)
colnames(grpid) <- c('Sample', 'MGrpid')
polori <- aggregate(polori_mc ~ Sample, data = dat_imp, mean)
colnames(polori) <- c('Sample', 'MPolori')
perc_stigma <- aggregate(perc_stigma_mc ~ Sample, data = dat_imp, mean)
colnames(perc_stigma) <- c('Sample', 'MPercStigma')
ins_stigma <- aggregate(ins_stigma ~ Sample, data = dat_imp, mean)
colnames(ins_stigma) <- c('Sample', 'MInsStigma')
gini <- aggregate(gini2016s ~ Sample, data = dat_imp, mean)
colnames(gini) <- c('Sample', 'Gini2016')
gdp <- aggregate(gdp2016s ~ Sample, data = dat_imp, mean)
colnames(gdp) <- c('Sample', 'GDPpc2016')

GenProp <- prop.table(table(dat_imp$Sample, dat_imp$gen), 1)
GenProp <- as.data.frame.matrix(GenProp)

colnames(GenProp) <- c('male', 'female', 'inter', 'trans', 'other')
GenProp$Sample <- row.names(GenProp)

cntrylvl_dat <- merge(wilac, colac, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, age, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, grpid, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, polori, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, GenProp, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, perc_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, ins_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gini, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gdp, by='Sample')

str(cntrylvl_dat)

country_priors_ins <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyMInsStigma2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyMInsStigma2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "MAge"),
              prior(normal(0, 10), class = "b", coef = "male"),
              prior(normal(0, 10), class = "b", coef = "female"),
              prior(normal(0, 10), class = "b", coef = "inter"),
              prior(normal(0, 10), class = "b", coef = "trans"),
              prior(normal(0, 10), class = "b", coef = "other"),
              prior(normal(0, 10), class = "b", coef = "MGrpid"),
              prior(normal(0, 10), class = "b", coef = "MPolori"),
              prior(normal(0, 10), class = "b", coef = "GDPpc2016"),
              prior(normal(0, 10), class = "b", coef = "Gini2016")
              )

country_priors_perc <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyMPercStigma2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyMPercStigma2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "MAge"),
              prior(normal(0, 10), class = "b", coef = "male"),
              prior(normal(0, 10), class = "b", coef = "female"),
              prior(normal(0, 10), class = "b", coef = "inter"),
              prior(normal(0, 10), class = "b", coef = "trans"),
              prior(normal(0, 10), class = "b", coef = "other"),
              prior(normal(0, 10), class = "b", coef = "MGrpid"),
              prior(normal(0, 10), class = "b", coef = "MPolori"),
              prior(normal(0, 10), class = "b", coef = "GDPpc2016"),
              prior(normal(0, 10), class = "b", coef = "Gini2016")
              )


## country level regression

wilac_ins_country <- brm(MWilac ~ poly(MInsStigma, 2, raw=FALSE) +
                             ## control variables
                             male +
                             female +
                             inter +
                             trans +
                             other +
                             MAge +
                             MGrpid +
                             MPolori +
                             Gini2016 +
                             GDPpc2016,
                         data = cntrylvl_dat,
                         save_pars = save_pars(all = TRUE),
                         prior = country_priors_ins)
summary(wilac_ins_country)


colac_ins_country <- brm(MColac ~ poly(MInsStigma, 2, raw=FALSE) +
                             ## control variables
                             male +
                             female +
                             inter +
                             trans +
                             other +
                             MAge +
                             MGrpid +
                             MPolori +
                             Gini2016 +
                             GDPpc2016,
                         data = cntrylvl_dat,
                         save_pars = save_pars(all = TRUE),
                         prior = country_priors_ins)
summary(colac_ins_country)


wilac_perc_country <- brm(MWilac ~ poly(MPercStigma, 2, raw=FALSE) +
                              ## control variables
                              male +
                              female +
                              inter +
                              trans +
                              other +
                              MAge +
                              MGrpid +
                              MPolori +
                              Gini2016 +
                              GDPpc2016,
                          data = cntrylvl_dat,
                          save_pars = save_pars(all = TRUE),
                          prior = country_priors_perc)
summary(wilac_perc_country)


colac_perc_country <- brm(MColac ~ poly(MPercStigma, 2, raw=FALSE) +
                              ## control variables
                              male +
                              female +
                              inter +
                              trans +
                              other +
                              MAge +
                              MGrpid +
                              MPolori +
                              Gini2016 +
                              GDPpc2016,
                          data = cntrylvl_dat,
                          save_pars = save_pars(all = TRUE),
                          prior = country_priors_perc)
summary(colac_perc_country)

tab_model(wilac_perc_country, colac_perc_country, show.se = TRUE)
tab_model(wilac_ins_country, colac_ins_country, show.se = TRUE)

##Linear models for comparison

country_priors_ins_linear <- c(prior(normal(0, 10), class = "Intercept"),
                               prior(normal(0, 10), class = "b", coef = "MInsStigma"),
                               prior(normal(0, 10), class = "b", coef = "MAge"),
                               prior(normal(0, 10), class = "b", coef = "male"),
                               prior(normal(0, 10), class = "b", coef = "female"),
                               prior(normal(0, 10), class = "b", coef = "inter"),
                               prior(normal(0, 10), class = "b", coef = "trans"),
                               prior(normal(0, 10), class = "b", coef = "other"),
                               prior(normal(0, 10), class = "b", coef = "MGrpid"),
                               prior(normal(0, 10), class = "b", coef = "MPolori"),
                               prior(normal(0, 10), class = "b", coef = "GDPpc2016"),
                               prior(normal(0, 10), class = "b", coef = "Gini2016")
                               )

country_priors_perc_linear <- c(prior(normal(0, 10), class = "Intercept"),
                                prior(normal(0, 10), class = "b", coef = "MPercStigma"),
                                prior(normal(0, 10), class = "b", coef = "MAge"),
                                prior(normal(0, 10), class = "b", coef = "male"),
                                prior(normal(0, 10), class = "b", coef = "female"),
                                prior(normal(0, 10), class = "b", coef = "inter"),
                                prior(normal(0, 10), class = "b", coef = "trans"),
                                prior(normal(0, 10), class = "b", coef = "other"),
                                prior(normal(0, 10), class = "b", coef = "MGrpid"),
                                prior(normal(0, 10), class = "b", coef = "MPolori"),
                                prior(normal(0, 10), class = "b", coef = "GDPpc2016"),
                                prior(normal(0, 10), class = "b", coef = "Gini2016")
                                )

wilac_ins_country_linear <- brm(MWilac ~ MInsStigma +
                             ## control variables
                             male +
                             female +
                             inter +
                             trans +
                             other +
                             MAge +
                             MGrpid +
                             MPolori +
                             Gini2016 +
                             GDPpc2016,
                             data = cntrylvl_dat,
                             save_pars = save_pars(all = TRUE),
                         prior = country_priors_ins_linear)

colac_ins_country_linear <- brm(MColac ~ MInsStigma +
                             ## control variables
                             male +
                             female +
                             inter +
                             trans +
                             other +
                             MAge +
                             MGrpid +
                             MPolori +
                             Gini2016 +
                             GDPpc2016,
                             data = cntrylvl_dat,
                             save_pars = save_pars(all = TRUE),
                         prior = country_priors_ins_linear)

wilac_perc_country_linear <- brm(MWilac ~ MPercStigma +
                             ## control variables
                             male +
                             female +
                             inter +
                             trans +
                             other +
                             MAge +
                             MGrpid +
                             MPolori +
                             Gini2016 +
                             GDPpc2016,
                             data = cntrylvl_dat,
                             save_pars = save_pars(all = TRUE),
                         prior = country_priors_perc_linear)

colac_perc_country_linear <- brm(MColac ~ MPercStigma +
                             ## control variables
                             male +
                             female +
                             inter +
                             trans +
                             other +
                             MAge +
                             MGrpid +
                             MPolori +
                             Gini2016 +
                             GDPpc2016,
                             data = cntrylvl_dat,
                             save_pars = save_pars(all = TRUE),
                         prior = country_priors_perc_linear)

bayes_factor(wilac_ins_country, wilac_ins_country_linear)
bayes_factor(colac_ins_country, colac_ins_country_linear)
bayes_factor(wilac_perc_country, wilac_perc_country_linear)
bayes_factor(colac_perc_country, colac_perc_country_linear)



##########
## Plots
##########

## Individual Levels

g1 <- plot(conditional_effects(wilac_perc, "perc_stigma_gmc"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Perceived Stigma") +
    ylab("Collective Actions Instentions")
g2 <- plot(conditional_effects(colac_perc, "perc_stigma_gmc"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Perceived Stigma") +
    ylab("Participation in Collective Actions")
ggsave("plots/perc_individual.png", arrangeGrob(g1, g2, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)

g3 <- plot(conditional_effects(wilac_ins, "ins_stigma"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Institutional Stigma") +
    ylab("Collective Actions Instentions")
g4 <- plot(conditional_effects(colac_ins, "ins_stigma"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Institutional Stigma") +
    ylab("Participation in Collective Actions")
ggsave("plots/ins_individual.png", arrangeGrob(g3, g4, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)

## Country Level

library(ggrepel)

g5 <- plot(conditional_effects(wilac_perc_country, "MPercStigma"), plot=FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    geom_count(aes(x=MPercStigma, y=MWilac),
               data=cntrylvl_dat, col='tomato3', alpha=.7, show.legend=F,
               inherit.aes=FALSE) +
    geom_text_repel(aes(x=MPercStigma, y=MWilac, label=Sample), max.overlaps=Inf,
                    data=cntrylvl_dat,
                    inherit.aes=FALSE) +
    xlab('Perceived Stigma') +
    ylab('Collective Action Intentions')

g6 <- plot(conditional_effects(colac_perc_country, "MPercStigma"), plot=FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    geom_count(aes(x=MPercStigma, y=MColac),
               data=cntrylvl_dat, col='tomato3', alpha=.7, show.legend=F,
               inherit.aes=FALSE) +
    geom_text_repel(aes(x=MPercStigma, y=MColac, label=Sample), max.overlaps=Inf,
                    data=cntrylvl_dat,
                    inherit.aes=FALSE) +
    xlab('Perceived Stigma') +
    ylab('Participation in Collective Action')


g7 <- plot(conditional_effects(wilac_ins_country, "MInsStigma"), plot=FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    geom_count(aes(x=MInsStigma, y=MWilac),
               data=cntrylvl_dat, col='tomato3', alpha=.7, show.legend=F,
               inherit.aes=FALSE) +
    geom_text_repel(aes(x=MInsStigma, y=MWilac, label=Sample), max.overlaps=Inf,
                    data=cntrylvl_dat,
                    inherit.aes=FALSE) +
    xlab('Institutional Stigma') +
    ylab('Collective Action Intentions')

g8 <- plot(conditional_effects(colac_ins_country, "MInsStigma"), plot=FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    geom_count(aes(x=MInsStigma, y=MColac),
               data=cntrylvl_dat, col='tomato3', alpha=.7, show.legend=F,
               inherit.aes=FALSE) +
    geom_text_repel(aes(x=MInsStigma, y=MColac, label=Sample), max.overlaps=Inf,
                    data=cntrylvl_dat,
                    inherit.aes=FALSE) +
    xlab('Institutional Stigma') +
    ylab('Participation in Collective Action')


ggsave("plots/perc_cntry.png", arrangeGrob(g5, g6, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)
ggsave("plots/ins_cntry.png", arrangeGrob(g7, g8, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)



#######################
## Mediation analysis
#######################


mx_priors <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "ins_stigma"),
              prior(normal(0, 10), class = "b", coef = "age_gmc"),
              prior(normal(0, 2), class = "b", coef = "gen2"),
              prior(normal(0, 2), class = "b", coef = "gen3"),
              prior(normal(0, 2), class = "b", coef = "gen4"),
              prior(normal(0, 2), class = "b", coef = "gen5"),
              prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
              prior(normal(0, 10), class = "b", coef = "polori_gmc"),
              prior(normal(0, 10), class = "b", coef = "gdp2016s"),
              prior(normal(0, 10), class = "b", coef = "gini2016s")
              )

ymx_priors <- c(prior(normal(0, 10), class = "Intercept"),
                prior(normal(0, 10), class = "b", coef = "polyperc_stigma_mc2rawEQFALSE1"),
                prior(normal(0, 10), class = "b", coef = "polyperc_stigma_mc2rawEQFALSE2"),
                prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE1"),
                prior(normal(0, 10), class = "b", coef = "polyins_stigma2rawEQFALSE2"),
                prior(normal(0, 10), class = "b", coef = "age_gmc"),
                prior(normal(0, 2), class = "b", coef = "gen2"),
                prior(normal(0, 2), class = "b", coef = "gen3"),
                prior(normal(0, 2), class = "b", coef = "gen4"),
                prior(normal(0, 2), class = "b", coef = "gen5"),
                prior(normal(0, 10), class = "b", coef = "grpid_gmc"),
                prior(normal(0, 10), class = "b", coef = "polori_gmc"),
                prior(normal(0, 10), class = "b", coef = "gdp2016s"),
                prior(normal(0, 10), class = "b", coef = "gini2016s")
                )

mxmodel <- brm_multiple(perc_stigma_mc ~ ins_stigma +
                            ## control variables
                            gen +
                            age_gmc +
                            grpid_gmc +
                            polori_gmc +
                            gini2016s +
                            gdp2016s +
                            (1 | Sample),
                        data = a.mids,
                        prior = mx_priors,
                        cores = 2,
                        file = 'models/mxmodel')
summary(mxmodel)

ymxmodel_co <- brm_multiple(colac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                                poly(ins_stigma, 2, raw = FALSE)  +
                                ## control variables
                                gen +
                                age_gmc +
                                grpid_gmc +
                                polori_gmc +
                                gini2016s +
                                gdp2016s +
                                (1 | Sample),
                            data = a.mids,
                            prior = ymx_priors,
                            cores = 2,
                            file = 'models/ymxmodel_co')
summary(ymxmodel_co)

ymxmodel_wi <-  brm_multiple(wilac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                                 poly(ins_stigma, 2, raw = FALSE)  +
                                 ## control variables
                                 gen +
                                 age_gmc +
                                 grpid_gmc +
                                 polori_gmc +
                                 gini2016s +
                                 gdp2016s +
                                 (1 | Sample),
                             data = a.mids,
                             prior = ymx_priors,
                             cores = 2,
                             file = 'models/ymxmodel_wi')
summary(ymxmodel_wi)


tab_model(mxmodel, show.se = TRUE)
tab_model(ymxmodel_wi, ymxmodel_co, show.se = TRUE)

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
a <- sort(a)
i1 <- sort(i1)
b1 <- sort(b1)
b2 <- sort(b2)

df <- data.frame()
means_co <- c()
xvals <- seq(-1.17764, 3.22305, length.out = 100)
for (val in xvals){
    b <- b1 + 2*b2*(i1 + a*val + sum(mean(dat$age_gmc, na.rm=T),
                                     mean(dat$grpid_gmc, na.rm=T),
                                     mean(dat$polori_gmc, na.rm=T),
                                     mean(dat$gini2016s, na.rm=T),
                                     mean(dat$gdp2016s, na.rm=T)))
    b <- sort(b)
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
                alpha=0.1,
                fill='blue') +
    geom_ribbon(aes(x=xvals, ymin=CIll, ymax=CIuu),
                alpha=0.1,
                fill='red') +
    xlab('Institutional Stigma') +
    ylab('Instantaneous Indirect Effect on Participation in collective actions') +
    ylim(-20, 20)

####
## Collective actions intentions
####

a <- as.vector(posterior_samples(mxmodel, 'ins_stigma', as.matrix = T))
i1 <- as.vector(posterior_samples(mxmodel, 'b_Intercept', as.matrix = T))
b1 <- as.vector(posterior_samples(ymxmodel_wi, 'polyperc_stigma_mc2rawEQFALSE1', as.matrix=T))
b2 <- as.vector(posterior_samples(ymxmodel_wi, 'polyperc_stigma_mc2rawEQFALSE2', as.matrix=T))
a <- sort(a)
i1 <- sort(i1)
b1 <- sort(b1)
b2 <- sort(b2)

df <- data.frame()
means_wi <- c()
xvals <- seq(-1.17764, 3.22305, length.out = 100)
for (val in xvals){
    b <- b1 + 2*b2*(i1 + a*val + sum(mean(dat$age_gmc, na.rm=T),
                                     mean(dat$grpid_gmc, na.rm=T),
                                     mean(dat$polori_gmc, na.rm=T),
                                     mean(dat$gini2016s, na.rm=T),
                                     mean(dat$gdp2016s, na.rm=T)))
    b <- sort(b)
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
                alpha=0.1,
                fill='blue') +
    geom_ribbon(aes(x=xvals, ymin=CIll, ymax=CIuu),
                alpha=0.1,
                fill='red') +
    xlab('Institutional Stigma') +
    ylab('Instantaneous Indirect Effect on Collective actions intentions') +
    ylim(-20, 20)

######################

ggsave("plots/medplot.png", arrangeGrob(medplot_wi, medplot_co, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)
