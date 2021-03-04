library(lme4)
library(ggplot2)
library(Amelia)
library(merTools)
library(sjPlot)
library(psych)
library(gridExtra)
library(mice)
library(miceadds)


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

## getting complete data

completedData <- complete(a.mids, 'long')

dat_imp <- aggregate(completedData[-c(1,2)], by=list(completedData$ID), FUN=mean)
dat_imp$Group.1 <- NULL
dat_imp$Sample <- NULL

dat_imp <- merge(dat[c('ID','Sample')], dat_imp, by='ID')

dat_imp$genF <- as.factor(dat_imp$gen)
dat_imp$gini2016S <- as.vector(scale(dat_imp$gini2016))
dat_imp$gdp2016S <- as.vector(scale(dat_imp$gdp2016))

mxmodel <- lmer(perc_stigma_mc ~ ins_stigma +
                    genF +
                    age_gmc +
                    grpid_gmc +
                    polori_gmc +
                    gini2016S +
                    gdp2016S +
                    (1 | Sample),
                data = dat_imp)
summary(mxmodel)

ymxmodel_co <- lmer(colac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                     poly(ins_stigma, 2, raw = FALSE)  +
                     as.factor(gen) +
                     age_gmc +
                     grpid_gmc +
                     polori_gmc +
                     scale(gini2016) +
                     scale(gdp2016) +
                     (1 | Sample),
                 data = dat_imp)
summary(ymxmodel_co)

ymxmodel_wi <- lmer(wilac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                     poly(ins_stigma, 2, raw = FALSE)  +
                     as.factor(gen) +
                     age_gmc +
                     grpid_gmc +
                     polori_gmc +
                     scale(gini2016) +
                     scale(gdp2016) +
                     (1 | Sample),
                 data = dat_imp)
summary(ymxmodel_wi)


## instantaneous indirect effects calculation

iie <- function(M, Y, x, data){
    M_est <- fixef(M)
    Y_est <- fixef(Y)
    control_M_est <- M_est[-c(1,2)]
    controlnames <- names(control_M_est)
    means <- rep(1, length=length(controlnames)) ## filling not found with 1s
    for (i in 1:length(controlnames)){
        try(
            means[i] <- mean(data[,controlnames[i]]),
            silent =  TRUE)
    }
    controls <- means*control_M_est
    a <- as.vector(M_est[2])
    b1 <- as.vector(Y_est[2])
    b2 <- as.vector(Y_est[3])
    i1 <- as.vector(M_est[1])
    theta <- a*(b1 + 2*b2*(i1 + a*x + sum(controls)))
    return(theta)
}


xvals <- seq(-1.17764, 3.22305, length.out = 100)
iievals_co <- iie(mxmodel, ymxmodel_co, xvals, dat_imp)
iievals_wi <- iie(mxmodel, ymxmodel_wi, xvals, dat_imp)

## Manual simple bootstrap


my_boot <- function(formulaM, formulaY, data, B=1000, xvar='ins_stigma', nx=100){

    n = dim(data)[1] 
    B = 1000
    xmin = min(data[xvar])
    xmax = max(data[xvar])

    results = matrix(nrow=B, ncol=nx)

    for(b in 1:B){
        i <- sample(x = 1:n, size = n, replace = TRUE)
        temp <- data[i,] 
        mxmodel_temp <- lmer(as.formula(formulaM),
                             data = temp)
        ymxmodel_temp <- lmer(as.formula(formulaY),
                              data = temp)
        xvals_temp <- seq(xmin, xmax, length.out = nx)
        iievals_temp <- iie(mxmodel_temp, ymxmodel_temp, xvals_temp, temp)
        results[b,] <- iievals_temp
    }
    results <- apply(results, 2, sort, decreasing=F)
    return(results)
}

fM <- perc_stigma_mc ~ ins_stigma +
                    genF +
                    age_gmc +
                    grpid_gmc +
                    polori_gmc +
                    gini2016S +
                    gdp2016S +
                    (1 | Sample)
fY_co <- colac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                     poly(ins_stigma, 2, raw = FALSE)  +
                     as.factor(gen) +
                     age_gmc +
                     grpid_gmc +
                     polori_gmc +
                     scale(gini2016) +
                     scale(gdp2016) +
                     (1 | Sample)
fY_wi <- wilac ~ poly(perc_stigma_mc, 2, raw = FALSE) +
                     poly(ins_stigma, 2, raw = FALSE)  +
                     as.factor(gen) +
                     age_gmc +
                     grpid_gmc +
                     polori_gmc +
                     scale(gini2016) +
                     scale(gdp2016) +
                     (1 | Sample)


results_co <- my_boot(fM, fY_co, dat_imp)
results_wi <- my_boot(fM, fY_wi, dat_imp)
results_wi

##plots colac

medplot_co <- ggplot() +
    geom_line(aes(x=xvals, y=iievals_co))+
    geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
    geom_ribbon(aes(x=xvals, y=iievals_co,
                    ymin=results_co[50,], ymax=results_co[950,]),
                alpha=0.1,
                fill='red')+
    geom_ribbon(aes(x=xvals, y=iievals_co,
                    ymin=results_co[25,], ymax=results_co[975,]),
                alpha=0.1,
                fill='blue')+
    geom_ribbon(aes(x=xvals, y=iievals_co,
                    ymin=results_co[10,], ymax=results_co[990,]),
                alpha=0.1,
                fill='black')+
    xlab('Institutional Stigma') +
    ylab('Instantaneous Indirect Effect on Participation in collective actions') +
    ylim(-20, 15)

##plots wilac

medplot_wi <- ggplot() +
    geom_line(aes(x=xvals, y=iievals_wi))+
    geom_abline(slope = 0, intercept = 0, linetype = "dashed")+
    geom_ribbon(aes(x=xvals, y=iievals_wi,
                    ymin=results_wi[50,], ymax=results_wi[950,]),
                alpha=0.1,
                fill='red')+
    geom_ribbon(aes(x=xvals, y=iievals_wi,
                    ymin=results_wi[25,], ymax=results_wi[975,]),
                alpha=0.1,
                fill='blue')+
    geom_ribbon(aes(x=xvals, y=iievals_wi,
                    ymin=results_wi[10,], ymax=results_wi[990,]),
                alpha=0.1,
                fill='black') +
    xlab('Institutional Stigma') +
    ylab('Instantaneous Indirect Effect on Collective actions intentions') +
    ylim(-20, 15)

ggsave("plots/medplot.png", arrangeGrob(medplot_co, medplot_wi, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)

## graphical explorations

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
    xlab("Perceived Stigma") +
    ylab("Collective Action Intentions")
g2 <- ggplot(dat_imp, aes(perc_stigma_gmc, res_colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Perceived Stigma") +
    ylab("Participation in Collective Actions")
ggsave("plots/perc_individual.png", arrangeGrob(g1, g2, nrow=1),width = 30, height = 12,
       units = "cm", limitsize = FALSE)

g3 <- ggplot(dat_imp, aes(ins_stigma, res_wilac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Institutional Stigma") +
    ylab("Collective Action Intentions")
g4 <- ggplot(dat_imp, aes(ins_stigma, res_colac)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    xlab("Institutional Stigma") +
    ylab("Participation in Collective Actions")
ggsave("plots/ins_individual.png", arrangeGrob(g3, g4, nrow=1),width = 30, height = 12,
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
summary(wilac_ins_country)


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
summary(colac_ins_country)


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
summary(wilac_perc_country)


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
summary(colac_perc_country)

tab_model(wilac_perc_country, colac_perc_country, show.ci = FALSE, show.se = TRUE)
tab_model(wilac_ins_country, colac_ins_country, show.ci = FALSE, show.se = TRUE)

library(ggrepel)

g5 <- ggplot(cntrylvl_dat, aes(MPercStigma, reswilac, label=Sample)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    geom_text_repel(aes(label=Sample), max.overlaps = Inf) +
    xlab('Perceived Stigma') +
    ylab('Collective Action Intentions') +
    ylim(-1.5, 1.5)

g6 <- ggplot(cntrylvl_dat, aes(MPercStigma, rescolac, label=Sample)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    geom_text_repel(aes(label=Sample), max.overlaps = Inf) +
    xlab('Perceived Stigma') +
    ylab('Participation in Collective Action') +
    ylim(-1.5, 1.5)

g7 <- ggplot(cntrylvl_dat, aes(MInsStigma, reswilac, label=Sample)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    geom_text_repel(aes(label=Sample), max.overlaps = Inf) +
    xlab('Institutional Stigma') +
    ylab('Collective Action Intentions') +
    ylim(-1.5, 1.5)

g8 <- ggplot(cntrylvl_dat, aes(MInsStigma, rescolac, label=Sample)) +
    geom_count(col="tomato3", show.legend=F) +
    stat_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
    geom_text_repel(aes(label=Sample), max.overlaps = Inf) +
    xlab('Institutional Stigma') +
    ylab('Participation in Collective Action') +
    ylim(-1.5, 1.5)

ggsave("plots/perc_cntry.png", arrangeGrob(g5, g6, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)
ggsave("plots/ins_cntry.png", arrangeGrob(g7, g8, nrow=1), width = 40, height = 20,
       units = "cm", limitsize = FALSE)
