install.packages('lme4')
install.packages('merTools')
install.packages('Amelia')
install.packages('lmerTest')
install.packages('mclust')

library(lme4)
library(ggplot2)
library(Amelia)
library(merTools)
library(lmerTest)
library(mice)
library(mclust)

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

dat_imp <- impute.out$imputations$imp1
dat_imp$ins_stigma <- as.numeric(dat_imp$ins_stigma)
str(dat_imp)


perc_stigma <- aggregate(perc_stigma_mc ~ Sample, data = dat_imp, mean)
ins_stigma <- aggregate(ins_stigma ~ Sample, data = dat_imp, mean)
age <- aggregate(age_mc ~ Sample, data = dat_imp, mean)
grpid <- aggregate(grpid_mc ~ Sample, data = dat_imp, mean)
gini2016 <- aggregate(gini2016 ~ Sample, data = dat_imp, mean)
gdp2016 <- aggregate(gdp2016 ~ Sample, data = dat_imp, mean)
wilac <- aggregate(wilac ~ Sample, data = dat_imp, mean)
colac <- aggregate(colac ~ Sample, data = dat_imp, mean)

countrydata <- merge(wilac, perc_stigma, by='Sample')
countrydata <- merge(ins_stigma, countrydata, by='Sample')
countrydata <- merge(colac, countrydata, by='Sample')
countrydata <- merge(age, countrydata, by='Sample')
countrydata <- merge(grpid, countrydata, by='Sample')
countrydata <- merge(gini2016, countrydata, by='Sample')
countrydata <- merge(gdp2016, countrydata, by='Sample')
countrydata$reswilac <- resid(lm(wilac ~  +
                                     age_mc +
                                     grpid_mc +
                                     scale(gdp2016) +
                                     scale(gini2016),
                                 data = countrydata))
countrydata$rescolac <- resid(lm(colac ~  +
                                     age_mc +
                                     grpid_mc +
                                     scale(gdp2016) +
                                     scale(gini2016),
                                 data = countrydata))



## country level regression

cntrywilac <- lm(wilac ~ poly(perc_stigma_mc, 2, raw=FALSE) +
                     age_mc +
                     grpid_mc +
                     scale(gdp2016) +
                     scale(gini2016),
                 data = countrydata)
summary(cntrywilac)

cntrycolac <- lm(colac ~ poly(perc_stigma_mc, 2, raw=FALSE) +
                     age_mc +
                     grpid_mc +
                     scale(gdp2016) +
                     scale(gini2016),
                 data = countrydata)
summary(cntrycolac)

cntrywilac <- lm(wilac ~ poly(ins_stigma, 2, raw=FALSE) +
                     age_mc +
                     grpid_mc +
                     scale(gdp2016) +
                     scale(gini2016),
                 data = countrydata)
summary(cntrywilac)

cntrycolac <- lm(colac ~ poly(ins_stigma, 2, raw=FALSE) +
                     age_mc +
                     grpid_mc +
                     scale(gdp2016) +
                     scale(gini2016),
                 data = countrydata)
summary(cntrycolac)


ggplot(countrydata, aes(perc_stigma_mc, reswilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

ggplot(countrydata, aes(perc_stigma_mc, rescolac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

ggplot(countrydata, aes(ins_stigma, reswilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1)

ggplot(countrydata, aes(ins_stigma, rescolac)) +
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
