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
library(mclust)
library(mediation)
library(sjPlot)
library(psych)
library(gridExtra)

dat <- read.csv('data/fulldataset.csv', row.names = 1)
names(dat)
str(dat)

dat$polori[dat$polori == 77] <- NA
dat$polori[dat$polori == 66] <- NA

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


## descriptives
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
dat$ins_stigma <- scale(dat$ins_stigma)

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

## multiple imputation step

set.seed(12345) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen'), idvars = c('Sample'), m = 5)
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
                              polori_gmc +
                              scale(gini2016) +
                              scale(gdp2016) +
                              (1 | Sample),
                          data = impute.out$imputations)
summary(wilac_perc)

##Willingness by institutional stigma model
wilac_ins <- lmerModList(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             as.factor(gen) +
                             age_gmc +
                             grpid_gmc +
                             polori_gmc +
                             scale(gini2016) +
                             scale(gdp2016) +
                             (1 | Sample),
                         data = impute.out$imputations)
summary(wilac_ins)

##Participation by perceive stigma model
colac_perc <- lmerModList(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                              ## control variables
                              as.factor(gen) +
                              age_gmc +
                              grpid_gmc +
                              polori_gmc +
                              scale(gini2016) +
                              scale(gdp2016) +
                              (1 | Sample),
                          data = impute.out$imputations)
summary(colac_perc)


##Participation by institutional stigma model
colac_ins <- lmerModList(colac ~ poly(ins_stigma, 2, raw=FALSE) +
                             ## control variables
                             as.factor(gen) +
                             age_gmc +
                             grpid_gmc +
                             polori_gmc +
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
                    polori_gmc +
                    scale(gini2016) +
                    scale(gdp2016) +
                    (1 | Sample),
                data = impute.out$imputations$imp1)

summary(mxmodel)

ymxmodel <- lmer(colac ~ perc_stigma_mc +
                     ins_stigma  +
                     as.factor(gen) +
                     age_gmc +
                     grpid_gmc +
                     polori_gmc +
                     scale(gini2016) +
                     scale(gdp2016) +
                     (1 | Sample),
                 data = impute.out$imputations$imp1)
summary(ymxmodel)

results <- mediation::mediate(mxmodel, ymxmodel, treat='ins_stigma', mediator='perc_stigma_mc',
                   boot=FALSE, sims = 100)
summary(results)

## graphical explorations

dat_imp <- impute.out$imputations$imp1

fit_explo_wilac <- lm(wilac ~ +
                    ## control variables
                    as.factor(gen) +
                    age_gmc +
                    grpid_gmc +
                    polori_gmc +
                    scale(gini2016) +
                    scale(gdp2016),
                    data = dat_imp)
fit_explo_colac <- lm(colac ~ +
                    ## control variables
                    as.factor(gen) +
                    age_gmc +
                    grpid_gmc +
                    polori_gmc +
                    scale(gini2016) +
                    scale(gdp2016),
                data = dat_imp)

dat_imp$res_wilac <- resid(fit_explo_wilac)
dat_imp$res_colac <- resid(fit_explo_colac)

g1 <- ggplot(dat_imp, aes(perc_stigma_gmc, res_wilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Perceive Stigma") +
    ylab("Collective Action Intentions")
g2 <- ggplot(dat_imp, aes(ins_stigma, res_wilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Institutional Stigma") +
    theme(axis.title.y = element_blank())
ggsave("plots/wilac_individual.png", arrangeGrob(g1, g2, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)

g3 <- ggplot(dat_imp, aes(perc_stigma_gmc, res_colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Perceive Stigma") +
    ylab("Participation in Collective Actions")
g4 <- ggplot(dat_imp, aes(ins_stigma, res_colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Institutional Stigma") +
    theme(axis.title.y = element_blank())
ggsave("plots/colac_individual.png", arrangeGrob(g3, g4, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)


## aggregated country level plots and analysis

wilac <- aggregate(wilac ~ Sample, data = dat_imp, mean)
colnames(wilac) <- c('Sample', 'MWilac')
colac <- aggregate(colac ~ Sample, data = dat_imp, mean)
colnames(colac) <- c('Sample', 'MColac')
reswilac <- aggregate(res_wilac ~ Sample, data = dat_imp, mean)
colnames(reswilac) <- c('Sample', 'reswilac')
rescolac <- aggregate(res_colac ~ Sample, data = dat_imp, mean)
colnames(rescolac) <- c('Sample', 'rescolac')
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
cntrylvl_dat <- merge(cntrylvl_dat, reswilac, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, rescolac, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, age, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, grpid, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, polori, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, GenProp, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, perc_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, ins_stigma, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gini, by='Sample')
cntrylvl_dat <- merge(cntrylvl_dat, gdp, by='Sample')

str(cntrylvl_dat)

## country level regression

options(browser="/usr/bin/google-chrome-stable")

wilac_ins_country <- lm(MWilac ~ poly(MInsStigma, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            MPolori +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(wilac_ins_country, use.viewer=FALSE)


colac_ins_country <- lm(MColac ~ poly(MInsStigma, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            MPolori +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(colac_ins_country, use.viewer=FALSE)


wilac_perc_country <- lm(MWilac ~ poly(MPercStigma, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            MPolori +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(wilac_perc_country, use.viewer=FALSE)


colac_perc_country <- lm(MColac ~ poly(MPercStigma, 2, raw=FALSE) +
                            ## control variables
                            male +
                            female+
                            inter+
                            trans+
                            other+
                            MAge +
                            MGrpid +
                            MPolori +
                            scale(Gini2016) +
                            scale(GDPpc2016),
                        data = cntrylvl_dat)
tab_model(colac_perc_country, use.viewer=FALSE)

tab_model(wilac_perc_country, colac_perc_country)
tab_model(wilac_ins_country, colac_ins_country)

g5 <- ggplot(cntrylvl_dat, aes(MPercStigma, reswilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab('Perceive Stigma') +
    ylab('Res. Collective Action Intentions') +
    ylim(-1.5, 1.5)

g6 <- ggplot(cntrylvl_dat, aes(MPercStigma, rescolac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab('Perceive Stigma') +
    ylab('Res. Participation in Collective Action') +
    ylim(-1.5, 1.5)

g7 <- ggplot(cntrylvl_dat, aes(MInsStigma, reswilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab('Institutional Stigma') +
    ylab('Res. Collective Action Intentions') +
    ylim(-1.5, 1.5)

g8 <- ggplot(cntrylvl_dat, aes(MInsStigma, rescolac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab('Institutional Stigma') +
    ylab('Res. Participation in Collective Action') +
    ylim(-1.5, 1.5)

ggsave("plots/perc_cntry.png", arrangeGrob(g5, g6, nrow=1), width = 30, height = 12,
       units = "cm", limitsize = FALSE)
ggsave("plots/ins_cntry.png", arrangeGrob(g7, g8, nrow=1), width = 30, height = 12,
       units = "cm", limitsize = FALSE)


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
