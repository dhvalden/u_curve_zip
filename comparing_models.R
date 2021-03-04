library(lme4)
library(ggplot2)
library(Amelia)
library(merTools)
library(sjPlot)
library(psych)
library(gridExtra)
library(mice)
library(miceadds)
library(lmtest)

options(browser="/usr/bin/google-chrome-stable")

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


## creating centered variables
dat <- get_mean_centered('perc_stigma', 'Sample', dat)
dat <- get_mean_centered('age', 'Sample', dat)
dat <- get_mean_centered('polori', 'Sample', dat)
dat <- get_mean_centered('grpid', 'Sample', dat)
dat$ins_stigma <- as.vector(scale(dat$ins_stigma))

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

set.seed(111) # set seed for reproducibility

impute.out <- amelia(dat, noms = c('gen'), idvars = c('Sample'), m = 5)
summary(impute.out)

a.mids <- datlist2mids(impute.out$imputations)

completedData <- complete(a.mids, 'long')

dat_imp <- aggregate(completedData[-c(1,2)], by=list(completedData$ID), FUN=mean)
dat_imp$Group.1 <- NULL
dat_imp$Sample <- NULL

dat_imp <- merge(dat[c('ID','Sample')], dat_imp, by='ID')


## individual level, random intercepts model

##Willingness by perceive stigma model
wilac_percP <- lmer(wilac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                        ## control variables
                        as.factor(gen) +
                        age_gmc +
                        grpid_gmc +
                        polori_gmc +
                        scale(gini2016) +
                        scale(gdp2016) +
                        (1 | Sample),
                    data = dat_imp)
wilac_percL <- lmer(wilac ~ perc_stigma_gmc +
                        ## control variables
                        as.factor(gen) +
                        age_gmc +
                        grpid_gmc +
                        polori_gmc +
                        scale(gini2016) +
                        scale(gdp2016) +
                        (1 | Sample),
                    data = dat_imp)

lrtest(wilac_percP, wilac_percL)

##Participation by perceive stigma model
colac_percP <- lmer(colac ~ poly(perc_stigma_gmc, 2, raw=FALSE) +
                        ## control variables
                        as.factor(gen) +
                        age_gmc +
                        grpid_gmc +
                        polori_gmc +
                        scale(gini2016) +
                        scale(gdp2016) +
                        (1 | Sample),
                    data = dat_imp)
colac_percL <- lmer(colac ~ perc_stigma_gmc +
                        ## control variables
                        as.factor(gen) +
                        age_gmc +
                        grpid_gmc +
                        polori_gmc +
                        scale(gini2016) +
                        scale(gdp2016) +
                        (1 | Sample),
                    data = dat_imp)
lrtest(colac_percP, colac_percL)

##Willingness by institutional stigma model
wilac_insP <- lmer(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       polori_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = dat_imp)
wilac_insL <- lmer(wilac ~ ins_stigma +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       polori_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = dat_imp)
lrtest(wilac_insP, wilac_insL)

##Participation by institutional stigma model
colac_insP <- lmer(colac ~ poly(ins_stigma, 2, raw=FALSE) +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       polori_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = dat_imp)
colac_insL <- lmer(colac ~ ins_stigma +
                       ## control variables
                       as.factor(gen) +
                       age_gmc +
                       grpid_gmc +
                       polori_gmc +
                       scale(gini2016) +
                       scale(gdp2016) +
                       (1 | Sample),
                   data = dat_imp)
lrtest(colac_insP, colac_insL)


## aggregated country level plots and analysis

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

## country level regression\

wilac_perc_countryP <- lm(MWilac ~ poly(MPercStigma, 2, raw=FALSE) +
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
wilac_perc_countryL <- lm(MWilac ~ MPercStigma +
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
lrtest(wilac_perc_countryP, wilac_perc_countryL)

colac_perc_countryP <- lm(MColac ~ poly(MPercStigma, 2, raw=FALSE) +
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
colac_perc_countryL <- lm(MColac ~ MPercStigma +
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
lrtest(colac_perc_countryP, colac_perc_countryL)


wilac_ins_countryP <- lm(MWilac ~ poly(MInsStigma, 2, raw=FALSE) +
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
wilac_ins_countryL <- lm(MWilac ~ MInsStigma +
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
lrtest(wilac_ins_countryP, wilac_ins_countryL)


colac_ins_countryP <- lm(MColac ~ poly(MInsStigma, 2, raw=FALSE) +
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
colac_ins_countryL <- lm(MColac ~ MInsStigma+
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
lrtest(colac_ins_countryP, colac_ins_countryL)
