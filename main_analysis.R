install.packages('lme4')
install.packages('merTools')
install.packages('Amelia')
install.packages('lmerTest')
install.packages('mclust')
install.packages('sjPlot')
install.packages('mediation')

library(lme4)
library(ggplot2)
library(Amelia)
library(merTools)
#library(lmerTest) # note: this packages interfere with mediation
library(mice)
library(mclust)
library(mediation)
library(sjPlot)
library(psych)

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
    return(tempdata)
}


## descriptives
options(browser="/usr/bin/firefox")

tab_df(aggregate(ID ~ Sample, data = dat, FUN = length))
prop.table(table(dat$gen))*100

frqdat <- dat
frqdat$gen  <- as.factor(frqdat$gen)
tab_df(describe(frqdat[-c(1,2)]))
tab_corr(frqdat[-c(1,2)])


## creating centered variables
dat <- get_mean_centered('perc_stigma', 'Sample', dat)
dat <- get_mean_centered('age', 'Sample', dat)
## dat <- get_mean_centered('polori', 'Sample', dat)
dat <- get_mean_centered('grpid', 'Sample', dat)
dat$ins_stigma <- scale(dat$ins_stigma)

str(dat)
names(dat)
summary(dat)

## checking correlations
round(cor(dat[,-c(1)], use = "pairwise"), 3)

## remove colliniar variables
dat$gai <- NULL
dat$gbgr <- NULL
dat$perc_stigma <- NULL
dat$age <- NULL
dat$grpid <- NULL
dat$perc_stigma_gm <- NULL
dat$age_gm <- NULL
dat$grpid_gm <- NULL

## multiple imputation step

impute.out <- amelia(dat, noms = c('gen'), idvars = c('Sample'), m = 20)
summary(impute.out)


## null models
summary(lmerModList(wilac ~ 1 + (1 | Sample),
                    data=impute.out$imputations))
summary(lmerModList(colac ~ 1 + (1 | Sample),
                    data=impute.out$imputations))

## individual level, random intercepts model

##Willingness by perceive stigma model
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

##Participation by perceive stigma model
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

##Willingness by institutional stigma model
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

##Participation by institutional stigma model
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

## Mediation analysis

mxmodel <- lmer(perc_stigma_mc ~ ins_stigma +
                    as.factor(gen) +
                    age_gmc +
                    grpid_gmc +
                    scale(gini2016) +
                    scale(gdp2016) +
                    (1 | Sample),
                data = impute.out$imputations$imp1)

summary(mxmodel)

ymxmodel <- lmer(wilac ~ perc_stigma_mc +
                     ins_stigma  +
                     as.factor(gen) +
                     age_gmc +
                     grpid_gmc +
                     scale(gini2016) +
                     scale(gdp2016) +
                     (1 | Sample),
                 data = impute.out$imputations$imp1)
summary(ymxmodel)

results <- mediate(mxmodel, ymxmodel, treat='ins_stigma', mediator='perc_stigma_mc',
                   boot=FALSE, sims = 100)
summary(results)

## graphical explorations

dat_imp <- impute.out$imputations$imp1

fit_explo_wilac <- lm(wilac ~ +
                    ## control variables
                    as.factor(gen) +
                    age_gmc +
                    grpid_gmc +
                    scale(gini2016) +
                    scale(gdp2016),
                    data = dat_imp)
fit_explo_colac <- lm(colac ~ +
                    ## control variables
                    as.factor(gen) +
                    age_gmc +
                    grpid_gmc +
                    scale(gini2016) +
                    scale(gdp2016),
                data = dat_imp)

dat_imp$res_wilac <- resid(fit_explo_wilac)
dat_imp$res_colac <- resid(fit_explo_colac)

ggplot(dat_imp, aes(perc_stigma_gmc, res_wilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Perceive Institutional Stigma") +
    ylab("Willingness to participate in activism") +
    ggtitle("Residualized polynomial regression plot")
ggsave("plots/wilac_perc.png", width = 20, height = 18,
       units = "cm", limitsize = FALSE)

ggplot(dat_imp, aes(perc_stigma_gmc, res_colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Perceive Institutional Stigma") +
    ylab("Participation in activism") +
    ggtitle("Residualized polynomial regression plot")
ggsave("plots/colac_perc.png", width = 20, height = 18,
       units = "cm", limitsize = FALSE)

ggplot(dat_imp, aes(ins_stigma, res_wilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Institutional Stigma") +
    ylab("Willingness to participate in activism") +
    ggtitle("Residualized polynomial regression plot")
ggsave("plots/wilac_ins.png", width = 20, height = 18,
       units = "cm", limitsize = FALSE)

ggplot(dat_imp, aes(ins_stigma, res_colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Institutional Stigma") +
    ylab("Participation in activism") +
    ggtitle("Residualized polynomial regression plot")
ggsave("plots/colac_ins.png", width = 20, height = 18,
       units = "cm", limitsize = FALSE)


## aggregated country level plots and analysis

dat_imp <- impute.out$imputations$imp10

wilac <- aggregate(wilac ~ Sample, data = dat_imp, mean)
colnames(wilac) <- c('Sample', 'MWilac')
colac <- aggregate(colac ~ Sample, data = dat_imp, mean)
colnames(colac) <- c('Sample', 'MColac')
age <- aggregate(age_mc ~ Sample, data = dat_imp, mean)
colnames(age) <- c('Sample', 'MAge')
grpid <- aggregate(grpid_mc ~ Sample, data = dat_imp, mean)
colnames(grpid) <- c('Sample', 'MGrpid')
perc_stigma <- aggregate(perc_stigma_mc ~ Sample, data = dat_imp, mean)
colnames(perc_stigma) <- c('Sample', 'MPercAccept')
ins_stigma <- aggregate(ins_stigma ~ Sample, data = dat_imp, mean)
colnames(ins_stigma) <- c('Sample', 'MInsAccept')
gini <- aggregate(gini2016 ~ Sample, data = dat_imp, mean)
colnames(gini) <- c('Sample', 'Gini2016')
gdp <- aggregate(gdp2016 ~ Sample, data = dat_imp, mean)
colnames(gdp) <- c('Sample', 'GDPpc2016')

GenProp <- prop.table(table(dat_imp$Sample, dat_imp$gen), 1)
GenProp <- as.data.frame.matrix(GenProp)
colnames(GenProp) <- c('male', 'female', 'inter', 'trans', 'other')
GenProp$Sample <- row.names(GenProp)
GenProp

cntrylvl_dat <- merge(wilac, colac, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, age, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, grpid, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, GenProp, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, perc_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, ins_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gini, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gdp, by='Sample')

str(cntrylvl_dat)

## country level regression

options(browser="/usr/bin/firefox")

wilac_ins_country <- lm(MWilac ~ poly(MInsAccept, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(wilac_ins_country, use.viewer=FALSE)


colac_ins_country <- lm(MColac ~ poly(MInsAccept, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(colac_ins_country, use.viewer=FALSE)



wilac_perc_country <- lm(MWilac ~ poly(MPercAccept, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(wilac_perc_country, use.viewer=FALSE)


colac_perc_country <- lm(MColac ~ poly(MPercAccept, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(colac_perc_country, use.viewer=FALSE)


ggplot(cntrylvl_dat, aes(MPercStigma, MWilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

ggplot(cntrylvl_dat, aes(perc_stigma_mc, rescolac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

ggplot(cntrylvl_dat, aes(ins_stigma, reswilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

ggplot(cntrylvl_dat, aes(ins_stigma, rescolac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) 





## models for clustering based on estimates of HLM

##Willingness by perceive stigma model
clust_wilac_perc <- lmer(wilac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              as.factor(gen) +
                              age_gmc +
                              grpid_gmc +
                              scale(gini2016) +
                              scale(gdp2016) +
                              (1 + poly(perc_stigma_gmc, 2, raw=FALSE) | Sample),
                          data = dat_imp)

coefs <- coef(clust_wilac_perc)
coefs

forclust <- as.data.frame(coefs$Sample['(Intercept)'])
colnames(forclust) <- c('intercept')
forclust$linearEffect <- coefs$Sample[,'poly(perc_stigma_gmc, 2, raw = FALSE)1']
forclust$quadraticEffect <- coefs$Sample[,'poly(perc_stigma_gmc, 2, raw = FALSE)2']
forclust$WillingessActivism_mean <- aggregate(wilac ~ Sample, mean, data=dat_imp)[,2]
forclust$InsAcceptance <- aggregate(ins_stigma ~ Sample, mean, data=dat_imp)[,2]
forclust$Gini2016 <- aggregate(gini2016 ~ Sample, mean, data=dat_imp)[,2]
forclust$GDPpc2016 <- aggregate(gdp2016 ~ Sample, mean, data=dat_imp)[,2]
forclust$PercAcceptance_cmean <- aggregate(perc_stigma_mc ~ Sample,
                                 mean, data=dat_imp)[,2]
forclust$sample <- row.names(forclust)
forclust
write.csv(forclust, 'data/wilac_clustering_data.csv')

##Participation by perceive stigma model
clust_colac_perc <- lmer(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              as.factor(gen) +
                              age_gmc +
                              grpid_gmc +
                              scale(gini2016) +
                              scale(gdp2016) +
                              (1 + poly(perc_stigma_gmc, 2, raw=FALSE) | Sample),
                          data = dat_imp)

coefs <- coef(clust_colac_perc)
coefs

forclust <- as.data.frame(coefs$Sample['(Intercept)'])
colnames(forclust) <- c('intercept')
forclust$linearEffect <- coefs$Sample[,'poly(perc_stigma_gmc, 2, raw = FALSE)1']
forclust$quadraticEffect <- coefs$Sample[,'poly(perc_stigma_gmc, 2, raw = FALSE)2']
forclust$ParticipationActivism_mean <- aggregate(colac ~ Sample, mean, data=dat_imp)[,2]
forclust$InsAcceptance <- aggregate(ins_stigma ~ Sample, mean, data=dat_imp)[,2]
forclust$Gini2016 <- aggregate(gini2016 ~ Sample, mean, data=dat_imp)[,2]
forclust$GDPpc2016 <- aggregate(gdp2016 ~ Sample, mean, data=dat_imp)[,2]
forclust$PercAcceptance_cmean <- aggregate(perc_stigma_mc ~ Sample,
                                 mean, data=dat_imp)[,2]
forclust$sample <- row.names(forclust)
forclust
write.csv(forclust, 'data/colac_clustering_data.csv')
