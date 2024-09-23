library(sp)
library(raster)
library(terra)

library(rjson)
library(geojsonR)

library(doMC)
library(doParallel)

###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(args[3])
# numSite <- 1


########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/mangrove/input/PLCM_Parameters.json')
source(params$setup$rFunctions)


########################################
## Get site name and image directory
geojsonDir <- params$setup$geojsonDir

strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckDir <- paste0(params$setup$outDir,strSite,'/chunk')
print(ckDir)



########################################
## 
imgBase <- raster(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- length(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks

nn <- 1
ptShp <- shapefile(paste0(params$setup$workDir,'shp/pts/mg_pts_',nn,'.shp'))
ptShp <- spTransform(ptShp,crs(imgBase))

pixNums <- extract(imgNum,ptShp)


##
dat <- matrix(NA,407,4)

for(ppp in 1:3){
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

  blue <- blue/10000; green <- green/10000; red <- red/10000; nir <- nir/10000
  
  pheno_pars <- params$phenology_parameters
  qa_pars    <- params$qa_parameters
  
  
  i1   <- (nir - red) / (nir + red) # NDVI
  i2   <- 2.5*(nir - red) / (nir + 2.4*red + 1) # EVI2
  i3   <- (nir - green) / (nir + green) # GNDVI
  i4   <- (green - nir) / (green + nir) # NDWI
  
  subdates <- dates[dates > 18071 & dates < 18627]
  subvi    <- i2[dates > 18071 & dates < 18627]
  
  if(ppp==1){
    dat[,1] <- subdates  
    dat[,2] <- subvi
  }else{
    dat[,ppp+1] <- subvi
  }
  
}


########################################
library(RColorBrewer)

mycol <- brewer.pal(3,)


setwd('/projectnb/modislc/users/mkmoon/mangrove/figures/')
png(filename=paste0('ts.png'),width=7.5,height=7.5,units='in',res=300)

par(mfrow=c(3,1),oma=c(2,1,0,0),mar=c(4,5,1,2),mgp=c(2.8,1.2,0))
plot(as.Date(dat[,1]),dat[,3],ylim=c(-0.05,0.75),pch=21,bg='red',cex=1.8,
     xlab='',ylab='EVI2',cex.lab=2,cex.axis=1.5)
plot(as.Date(dat[,1]),dat[,4],ylim=c(-0.05,0.75),pch=21,bg='darkgreen',cex=1.8,
     xlab='',ylab='EVI2',cex.lab=2,cex.axis=1.5)
plot(as.Date(dat[,1]),dat[,2],ylim=c(-0.05,0.75),pch=21,bg='darkorange',cex=1.8,
     xlab='',ylab='EVI2',cex.lab=2,cex.axis=1.5)

dev.off()



