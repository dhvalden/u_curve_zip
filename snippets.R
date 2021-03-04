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
