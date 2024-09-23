library(sp)
library(raster)
library(terra)
library(sf)

library(rjson)
library(geojsonR)

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
ppp <- 1
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


########################################
GetFeatures <- function(blue, green, red, nir, dates, phenYrs, params){
  
  # Despike, calculate dormant value, fill negative VI values with dormant value
  log <- try({    
    
    pheno_pars <- params$phenology_parameters
    qa_pars    <- params$qa_parameters
    
    b2 <- blue/10000; b3 <- green/10000; b4 <- red/10000; b5 <- nir/10000
    
    i1   <- (b5 - b4) / (b5 + b4) # NDVI
    i2   <- (b5 - b3) / (b5 + b3) # GNDVI
    i3   <- b5 * i1               # NIRv
    vi   <- 2.5*(b5 - b4) / (b5 + 2.4*b4 + 1) # EVI2
    
    
    # Spikes check, and remove
    spikes     <- CheckSpike_MultiBand(b2, b4, vi, dates, pheno_pars)
    b2[spikes] <- NA; b3[spikes] <- NA; b4[spikes] <- NA; b5[spikes] <- NA
    i1[spikes] <- NA; i2[spikes] <- NA; i3[spikes] <- NA; vi[spikes] <- NA
    
    # Replace negative VIs with dormant value
    dormIms <- dates >= pheno_pars$dormStart & dates <= pheno_pars$dormEnd
    vi_dorm <- quantile(vi[dormIms & vi>0],probs=pheno_pars$dormantQuantile,na.rm=T)   # Calc vi dormant value using non-negative VIs
    
    
    #now calculate dormancy values and fill individual bands
    dormObs <- dormIms & vi < vi_dorm    #Defining dormant observations for bands as median on dates when vi < vi_dorm
    b2_dorm <- median(b2[dormObs], na.rm=T); b2[dormObs] <- b2_dorm
    b3_dorm <- median(b3[dormObs], na.rm=T); b3[dormObs] <- b3_dorm
    b4_dorm <- median(b4[dormObs], na.rm=T); b4[dormObs] <- b4_dorm
    b5_dorm <- median(b5[dormObs], na.rm=T); b5[dormObs] <- b5_dorm
    i1_dorm <- median(i1[dormObs], na.rm=T); i1[dormObs] <- i1_dorm
    i2_dorm <- median(i1[dormObs], na.rm=T); i2[dormObs] <- i2_dorm
    i3_dorm <- median(i1[dormObs], na.rm=T); i3[dormObs] <- i3_dorm
    vi[vi < vi_dorm] <- vi_dorm
    
    
    #
    splineStart <- as.Date(as.Date(paste0(phenYrs,'-01-01')) - pheno_pars$splineBuffer) 
    numDaysFit  <- 365 + (pheno_pars$splineBuffer * 2)    
    splineEnd   <- splineStart+(numDaysFit-1)
    all_dates   <- seq(min(splineStart), max(splineEnd), by="day")
    
    daysVec <- 1:numDaysFit
    inYear  <- daysVec > pheno_pars$splineBuffer & daysVec <= (pheno_pars$splineBuffer+365)
    prevYear <- daysVec <= pheno_pars$splineBuffer
    nextYear <- daysVec > (pheno_pars$splineBuffer+365)
    
    numYrs <- length(phenYrs)
    vecLength <- numDaysFit*numYrs
    
    #
    smoothMat <- matrix(NA, numDaysFit, numYrs)
    b2Mat <- smoothMat; b3Mat <- smoothMat; b4Mat <- smoothMat; b5Mat <- smoothMat
    i1Mat <- smoothMat; i2Mat <- smoothMat; i3Mat <- smoothMat
    
  },silent=TRUE)
  #If there is an error despiking or other initial steps, return NAs
  if(inherits(log, "try-error")){return(matrix(NA,58*length(phenYrs)))}   
  
    
  outAll <- c()
  for(y in 1:numYrs){
    log <- try({    
      
      pred_dates <- seq(splineStart[y], splineEnd[y], by="day")
      
      dateRange <- dates >= splineStart[y] & dates <= splineEnd[y] & !is.na(vi)   
      dateSub <- dates[dateRange]
      viSub <- vi[dateRange]
      
      smoothed_vi <- Smooth_VI(viSub, dateSub, pred_dates, pheno_pars, vi_dorm)  #Fit spline
      
      # Number of clear observation
      numObs <- sum(dateSub > splineStart[y]+184 & dateSub < splineEnd[y]-184)   #Number of observations in year
      
      
      # Fit spline to individual bands
      smooth_b2 <- Smooth_Bands(b2[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_b3 <- Smooth_Bands(b3[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_b4 <- Smooth_Bands(b4[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_b5 <- Smooth_Bands(b5[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_i1 <- Smooth_Bands(i1[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_i2 <- Smooth_Bands(i2[dateRange], dateSub, pred_dates, pheno_pars)
      smooth_i3 <- Smooth_Bands(i3[dateRange], dateSub, pred_dates, pheno_pars)
      
      smoothed_bns <- cbind(smooth_b2, smooth_b3, smooth_b4, smooth_b5, smooth_i1, smooth_i2, smooth_i3, smoothed_vi)
      
      
      
      ################################################
      #Fit phenology
      peaks <- FindPeaks(smoothed_vi)
      if (all(is.na(peaks))) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #Find full segments
      full_segs <- GetSegs(peaks, smoothed_vi, pheno_pars)
      if (is.null(full_segs)) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #Only keep segments with peaks within year *****
      full_segs <- full_segs[inYear[sapply(full_segs, "[[", 2)] ]  #check if peaks are in the year
      if (length(full_segs)==0) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #Get PhenoDates
      pheno_dates <- GetPhenoDates(full_segs, smoothed_vi, pred_dates, pheno_pars)
      phen <- unlist(pheno_dates, use.names=F)
      phen <- phen - as.numeric(as.Date(paste0((as.numeric(phenYrs[y])-1),'-12-31')))
      if (all(is.na(phen))) {outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      #EVI layers
      seg_metrics <- lapply(full_segs, GetSegMetrics,smoothed_vi,viSub,pred_dates,dateSub) #full segment metrics
      un <- unlist(seg_metrics, use.names=F)
      ln <- length(un)
      seg_amp <- un[seq(1, ln, by=9)] * 10000
      seg_max <- un[seq(2, ln, by=9)] * 10000
      seg_int <- un[seq(3, ln, by=9)] * 100  
      gup_rsq <- un[seq(4, ln, by=9)] * 10000
      gup_maxgap <- un[seq(6, ln, by=9)] 
      gdown_rsq <- un[seq(7, ln, by=9)] * 10000
      gdown_maxgap <- un[seq(9, ln, by=9)]
      
      
      ##
      theOrd <- order(seg_amp,decreasing=T)   
      
      # Filter for bad EVI layers
      if(seg_max[theOrd[1]] > 10000 | seg_max[theOrd[1]] < 0 | seg_amp[theOrd[1]] > 10000 | seg_amp[theOrd[1]] < 0){
        outAll <- c(outAll,rep(NA,58));next}
      
      numRecords <- length(seg_amp)  #how many cycles were recorded
      naCheck <- is.na(seg_amp)
      numCyc <- sum(naCheck == 0)  #how many cycles have good data (seg metrics has valid observations)
      
      
      # QA
      if(length(full_segs)==1){
        qual_1 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)
      }else{
        qual_1 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)[[1]][1]
        qual_2 <- GetQAs(gup_rsq, gdown_rsq, gup_maxgap, gdown_maxgap, theOrd, qa_pars)[[2]][1]
      } 
      
      
      ################################################
      if(numCyc == 0){outAll <- c(outAll,annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),rep(NA,18));next}
      
      if(numRecords == 1) {
        out <- c(annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),
                 1,phen,qual_1,c(rep(NA,7),4),numObs)
      }else{
        phen1 <- phen[seq(theOrd[1], length(phen), by = numRecords)]
        phen2 <- phen[seq(theOrd[2], length(phen), by = numRecords)]
        if(naCheck[theOrd[2]]){
          out <- c(annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),
                   numCyc,phen1,qual_1,c(rep(NA,7),4),numObs)  
        }else{
          out <- c(annualMetrics(viSub,dateSub,smoothed_bns,pred_dates,phenYrs[y],pheno_pars,vi_dorm),
                   numCyc,phen1,qual_1,phen2,qual_2,numObs)  
        }
      }
      
    },silent=TRUE) #End of the try block
    if(inherits(log, "try-error")){outAll <- c(outAll,rep(NA,58))
    }else{outAll <- c(outAll,out);remove(out)}
    
  }
  
  return(outAll)
}



