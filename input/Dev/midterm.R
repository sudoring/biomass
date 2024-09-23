library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 1


###############################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/mangrove/input/PLCM_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name, image directory and coordinate
geojsonDir <- params$setup$geojsonDir

siteInfo <- GetSiteInfo(numSite,geojsonDir,params)

imgDir <- siteInfo[[1]]; strSite <- siteInfo[[2]]
print(paste(strSite,';',imgDir))

cLong <- siteInfo[[3]];cLat <- siteInfo[[4]]
print(paste(cLong,';',cLat))
cLong <- 115.20700
cLat  <- -8.751338


########################################
### Aqusition dates
## List of files
dfileSR   <- list.files(path=imgDir,pattern=glob2rx('*MS_SR*.tif'),recursive=T)
dfileUDM2 <- list.files(path=imgDir,pattern=glob2rx('*_udm2*.tif'),recursive=T)

fileSR   <- list.files(path=imgDir,pattern=glob2rx('*MS_SR*.tif'),recursive=T,full.names=T)
fileUDM2 <- list.files(path=imgDir,pattern=glob2rx('*_udm2*.tif'),recursive=T,full.names=T)


# Dates
yy <- c(); mm <- c(); dd <- c()
for(i in 1:length(dfileSR)){
  yy[i] <- substr(unlist(strsplit(dfileSR[i],'/'))[4],3,4)
  mm[i] <- substr(unlist(strsplit(dfileSR[i],'/'))[4],5,6)
  dd[i] <- substr(unlist(strsplit(dfileSR[i],'/'))[4],7,8)
  
  dfileSR[i] <- unlist(strsplit(dfileSR[i],'/'))[4]
  dfileUDM2[i] <- unlist(strsplit(dfileUDM2[i],'/'))[4]
}
dates_all <- as.Date(paste(mm,'/',dd,'/',yy,sep=''),'%m/%d/%y')
dates <- unique(dates_all)

#
datesod <- order(dates)
setwd('/projectnb/modislc/users/mkmoon/mangrove/figures/')
png(filename='numofimage.png',width=7.5,height=6.5,unit='in',res=300)
par(oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
plot(dates[datesod],
     xlab='Number of Image',ylab='Dates',cex.axis=1.2,cex.lab=1.5)
dev.off()
#

print(length(dates))






########################################
numSite <- 1


########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/mangrove/input/PLCM_Parameters.json')
source(params$setup$rFunctions)

## Get site name and image directory
geojsonDir <- params$setup$geojsonDir

strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)


########################################
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks

##
ptShp <- shapefile(paste0(params$setup$workDir,'shp/mg_pts_',numSite,'.shp'))
ptShp <- spTransform(ptShp,crs(imgBase))
pixNums <- extract(imgNum,ptShp)

##
ppp <- 4
pixNum  <- pixNums[ppp]

ckNum <- sprintf('%03d',(pixNum%/%chunk+1))
file <- list.files(path=ckDir,pattern=glob2rx(paste0('*',ckNum,'.rda')),full.names=T)

load(file)

##
blue  <- band1[pixNum%%chunk,]
green <- band2[pixNum%%chunk,]
red   <- band3[pixNum%%chunk,]
nir   <- band4[pixNum%%chunk,]
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr


b2 <- blue/10000; b3 <- green/10000; b4 <- red/10000; b5 <- nir/10000

i1   <- (b5 - b4) / (b5 + b4) # NDVI
i2   <- (b5 - b3) / (b5 + b3) # GNDVI
i3   <- b5 * i1               # NIRv
vi   <- 2.5*(b5 - b4) / (b5 + 2.4*b4 + 1) # EVI2


###
library(RColorBrewer)
mycol <- brewer.pal(8,'Set1')


setwd('/projectnb/modislc/users/mkmoon/mangrove/figures/')
png(filename='ts_bns_3.png',width=11,height=6.5,unit='in',res=300)

par(oma=c(1,1,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(dates,vi,ylim=c(-0.1,1),pch=21,bg=mycol[1],
     xlab='Date',ylab='Indicies',cex=1.2,cex.lab=1.5,cex.axis=1.5)
points(dates,i1,pch=21,bg=mycol[2],cex=1.2)
points(dates,i2,pch=21,bg=mycol[3],cex=1.2)
points(dates,i3,pch=21,bg=mycol[4],cex=1.2)
points(dates,b2,pch=21,bg=mycol[5],cex=1.2)
points(dates,b3,pch=21,bg=mycol[6],cex=1.2)
points(dates,b4,pch=21,bg=mycol[7],cex=1.2)
points(dates,b5,pch=21,bg=mycol[8],cex=1.2)

legend('topleft',c('EVI2','NDVI','GNDVI','NIRv','Blue','Green','Red','NIR'),
       pch=21,pt.bg=mycol,cex=1.3)
dev.off()
