rm(list = ls())

library(terra)
library(sf)


############################################################
plt <- rast('/projectnb/modislc/users/mkmoon/biomass/data/taehwa_psscene_analytic_sr_udm2/PSScene/20240318_012943_54_24c0_3B_AnalyticMS_SR_clip.tif')

plt_r <- plt[[3]]/10000
plt_n <- plt[[4]]/10000 
evi2  <- 2.5 * (plt_n - plt_r) / (plt_n + 2.4 * plt_r + 1)

plot(evi2)


# gldRep <- project(gld,plt)
# gldRep <- crop(gldRep,plt)
# 
# plot(gldRep)


############################################################
outDir <- '/projectnb/modislc/users/mkmoon/biomass/data/'
writeRaster(evi2,filename=paste0(outDir,'evi2_1.tif'),overwrite=TRUE)


############################################################
plt <- rast('/projectnb/modislc/users/mkmoon/biomass/data/taehwa_1_psscene_analytic_sr_udm2/PSScene/20240418_021409_08_241c_3B_AnalyticMS_SR_clip.tif')

plt_r <- plt[[3]]/10000
plt_n <- plt[[4]]/10000 
evi2  <- 2.5 * (plt_n - plt_r) / (plt_n + 2.4 * plt_r + 1)

plot(evi2)

outDir <- '/projectnb/modislc/users/mkmoon/biomass/data/'
writeRaster(evi2,filename=paste0(outDir,'evi2_2.tif'),overwrite=TRUE)
