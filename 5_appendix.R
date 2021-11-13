library(ggplot2)
library(Amelia)
library(sjPlot)
library(psych)
library(gridExtra)
library(mice)
library(miceadds)
library(brms)

 options(browser="firefox")

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
             'wilac_l', ## low cost
             'colac_l', ## low cost
             'ins_stigma',
             'perc_stigma',
             'gini2016',
             'gdp2016')]

str(dat)

## renaming low cost collective actions for convenince
colnames(dat)[colnames(dat) == 'wilac_l'] <- 'wilac'
colnames(dat)[colnames(dat) == 'colac_l'] <- 'colac'

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
dat <- get_mean_centered('polori', 'Sample', dat)
dat <- get_mean_centered('grpid', 'Sample', dat)
dat$ins_stigma <- as.vector(scale(dat$ins_stigma))
dat$gini2016s <- as.vector(scale(dat$gini2016))
dat$gdp2016s <- as.vector(scale(dat$gdp2016))


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
dat$new_gen <- as.factor(dat$new_gen)
set.seed(111) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen_min', 'new_gen'), idvars = c('Sample'), m = 5)
summary(impute.out)
 
a.mids <- datlist2mids(impute.out$imputations)

########################################
## Main Bayesian regression analysis
########################################

####
## Individual Level
####

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
                           file = 'models/wilac_perc_low'
                           )
summary(wilac_perc)

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
                          file = 'models/colac_perc_low'
                          )
summary(colac_perc)

##Willingness by institutional stigma model
wilac_ins <- brm_multiple(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
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
                         prior = mypriors_ins,
                         cores = 3,
                         seed = 123,
                         refresh = 0,
                         open_progress = FALSE,
                         save_pars = save_pars(all = TRUE),
                         file = 'models/wilac_ins_low'
                         )
summary(wilac_ins)

##Participation by institutional stigma model

colac_ins <- brm_multiple(colac ~ poly(ins_stigma, 2, raw=FALSE) +
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
                         prior = mypriors_ins,
                         cores = 3,
                         seed = 123,
                         refresh = 0,
                         open_progress = FALSE,
                         save_pars = save_pars(all = TRUE),
                         file = 'models/colac_ins_low'
                         )
summary(colac_ins)

tab_model(wilac_perc, colac_perc, robust = TRUE)
tab_model(wilac_ins, colac_ins, robust = TRUE)

bayestestR::hdi(wilac_perc)
bayestestR::hdi(colac_perc)
bayestestR::hdi(wilac_ins)
bayestestR::hdi(colac_ins)

####
##Country Level
####

## getting complete data

completedData <- complete(a.mids, 'long')
dat_imp <- aggregate(completedData[-c(1,2)], by=list(completedData$ID), FUN=mean)
dat_imp$Group.1 <- NULL
dat_imp$Sample <- NULL
dat_imp$gen_min <- NULL
dat_imp$new_gen <- NULL

dat_imp <- merge(dat[c('ID','Sample', 'gen_min', 'new_gen')], dat_imp, by='ID')

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

Gen_minProp <- prop.table(table(dat_imp$Sample, dat_imp$gen_min), 1)
Gen_minProp <- as.data.frame.matrix(Gen_minProp)
colnames(Gen_minProp) <- c('prop_no_gen_min', 'prop_gen_min')
Gen_minProp$Sample <- row.names(Gen_minProp)

new_gen_prop <- prop.table(table(dat_imp$Sample, dat_imp$new_gen), 1)
new_gen_prop <- as.data.frame.matrix(new_gen_prop)
colnames(new_gen_prop) <- c('prop_male', 'prop_female', 'prop_other')
new_gen_prop$Sample <- row.names(new_gen_prop)

cntrylvl_dat <- merge(wilac, colac, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, age, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, grpid, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, polori, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, Gen_minProp, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, new_gen_prop, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, perc_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, ins_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gini, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gdp, by='Sample')

country_priors_ins <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyMInsStigma2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyMInsStigma2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "MAge"),
              prior(normal(0, 10), class = "b", coef = "prop_gen_min"),
              prior(normal(0, 10), class = "b", coef = "prop_female"),
              prior(normal(0, 10), class = "b", coef = "prop_other"),
              prior(normal(0, 10), class = "b", coef = "MGrpid"),
              prior(normal(0, 10), class = "b", coef = "MPolori"),
              prior(normal(0, 10), class = "b", coef = "GDPpc2016"),
              prior(normal(0, 10), class = "b", coef = "Gini2016")
              )

country_priors_perc <- c(prior(normal(0, 10), class = "Intercept"),
              prior(normal(0, 10), class = "b", coef = "polyMPercStigma2rawEQFALSE1"),
              prior(normal(0, 10), class = "b", coef = "polyMPercStigma2rawEQFALSE2"),
              prior(normal(0, 10), class = "b", coef = "MAge"),
              prior(normal(0, 10), class = "b", coef = "prop_gen_min"),
              prior(normal(0, 10), class = "b", coef = "prop_female"),
              prior(normal(0, 10), class = "b", coef = "prop_other"),
              prior(normal(0, 10), class = "b", coef = "MGrpid"),
              prior(normal(0, 10), class = "b", coef = "MPolori"),
              prior(normal(0, 10), class = "b", coef = "GDPpc2016"),
              prior(normal(0, 10), class = "b", coef = "Gini2016")
              )


## country level regression

wilac_ins_country <- brm(MWilac ~ poly(MInsStigma, 2, raw=FALSE) +
                             ## control variables
                             prop_gen_min +
                             prop_female +
                             prop_other +
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
                             prop_gen_min +
                             prop_female +
                             prop_other +
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
                             prop_gen_min +
                             prop_female +
                             prop_other +
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
                             prop_gen_min +
                             prop_female +
                             prop_other +
                             MAge +
                             MGrpid +
                             MPolori +
                             Gini2016 +
                             GDPpc2016,
                          data = cntrylvl_dat,
                          save_pars = save_pars(all = TRUE),
                          prior = country_priors_perc)
summary(colac_perc_country)

tab_model(wilac_perc_country, colac_perc_country, robust = TRUE)
tab_model(wilac_ins_country, colac_ins_country, robust = TRUE)

bayestestR::hdi(wilac_perc_country)
bayestestR::hdi(colac_perc_country)
bayestestR::hdi(wilac_ins_country)
bayestestR::hdi(colac_ins_country)

##########
## Plots
##########

## Individual Levels

g1 <- plot(conditional_effects(wilac_perc, "perc_stigma_gmc"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Perceived Stigma") +
    ylab("Low Cost Collective Actions Intentions")
g2 <- plot(conditional_effects(colac_perc, "perc_stigma_gmc"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Perceived Stigma") +
    ylab("Participation in Low Cost Collective Actions")
ggsave("plots/perc_individual_low.png", arrangeGrob(g1, g2, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)

g3 <- plot(conditional_effects(wilac_ins, "ins_stigma"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Institutional Stigma") +
    ylab("Low Cost Collective Actions Intentions")
g4 <- plot(conditional_effects(colac_ins, "ins_stigma"), plot = FALSE)[[1]] +
    theme_bw() +
    ylim(1, 7) +
    xlab("Institutional Stigma") +
    ylab("Participation in Low Cost Collective Actions")
ggsave("plots/ins_individual_low.png", arrangeGrob(g3, g4, nrow=1),width = 30, height = 12,
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
    ylab('Low Cost Collective Action Intentions')

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
    ylab('Participation in Low Cost Collective Action')


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
    ylab('Low Cost Collective Action Intentions')

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
    ylab('Participation in Low Cost Collective Action')


ggsave("plots/perc_cntry_low.png", arrangeGrob(g5, g6, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)
ggsave("plots/ins_cntry_low.png", arrangeGrob(g7, g8, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)
