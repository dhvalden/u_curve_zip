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

#############################################
## Within models low cost collective actions
#############################################

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


tab_model(wilac_perc, colac_perc, robust = TRUE)

bayestestR::hdi(wilac_perc)
bayestestR::hdi(colac_perc)

##############################################
## Between Models low cost collective actions
##############################################

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
                             file = 'models/wilac_perc_g_low'
                             )
summary(wilac_perc_g)

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
                             file = 'models/colac_perc_g_low'
                             )
summary(colac_perc_g)

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
                            file = 'models/wilac_ins_g_low'
                            )
summary(wilac_ins_g)


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
                            file = 'models/colac_ins_g_low'
                            )
summary(colac_ins_g)

bayestestR::hdi(wilac_perc_g)
bayestestR::hdi(colac_perc_g)
bayestestR::hdi(wilac_ins_g)
bayestestR::hdi(colac_ins_g)


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
