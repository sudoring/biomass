library(sp)
library(raster)
library(terra)
library(sf)
library(rjson)

library(randomForest)



###############################
args <- commandArgs()
print(args)

numSite <- as.numeric(substr(args[3],1,3))
yToDo   <- as.numeric(substr(args[3],4,7))
# numSite <- 1; yToDo <- 4



########################################
params <- fromJSON(file='/usr3/graduate/mkmoon/GitHub/mangrove/input/PLCM_Parameters.json')
source(params$setup$rFunctions)

productTable <- read.csv(params$setup$productTable,header=T,stringsAsFactors = F)
phenYrs <- params$setup$phenStartYr:params$setup$phenEndYr



########################################
## Get site name, image directory and coordinate
strSite <- list.dirs(params$setup$outDir,full.names=F,recursive=F)[numSite]
print(strSite)

ckPheDir <- paste0(params$setup$outDir,strSite,'/chunk_feat')
print(ckPheDir)

files <- list.files(path=ckPheDir,pattern=glob2rx('*.rda'),full.names=T)
print(length(files))



########################################
imgBase <- rast(paste0(params$setup$outDir,strSite,'/base_image.tif'))
numPix <- ncol(imgBase)*nrow(imgBase)
imgNum <- setValues(imgBase,1:numPix)
numChunks <- params$setup$numChunks
chunk <- numPix%/%numChunks



########################################
# Load training data
t1 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_1/training.shp')   # 1 mangrove
t2 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_2/training_2.shp') # 2 deep water
t3 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_3/training_3.shp') # 3 shallow water
t4 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_4/training_4.shp') # 4 terrestrial vegetation
t5 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_5/training_5.shp') # 5 buildings
t6 <- vect('/projectnb/modislc/users/mkmoon/mangrove/PLCM/shp/training_6/training_6.shp') # 6 crop

t1 <- project(t1,crs(imgBase)); t2 <- project(t2,crs(imgBase)); t3 <- project(t3,crs(imgBase))
t4 <- project(t4,crs(imgBase)); t5 <- project(t5,crs(imgBase)); t6 <- project(t6,crs(imgBase))

pt1 <- extract(imgNum,t1); pt2 <- extract(imgNum,t2); pt3 <- extract(imgNum,t3)
pt4 <- extract(imgNum,t4); pt5 <- extract(imgNum,t5); pt6 <- extract(imgNum,t6)


# Divide the points into training, validation, and testing; 8:1:1
set.seed(456789)
pt11s <- sample(1:nrow(pt1),round(nrow(pt1)*0.8)); pt11  <- pt1[pt11s,]; pt12   <- pt1[-pt11s,] # T1
pt21s <- sample(1:nrow(pt2),round(nrow(pt2)*0.8)); pt21  <- pt2[pt21s,]; pt22   <- pt2[-pt21s,] # T2
pt31s <- sample(1:nrow(pt3),round(nrow(pt3)*0.8)); pt31  <- pt3[pt31s,]; pt32   <- pt3[-pt31s,] # T3
pt41s <- sample(1:nrow(pt4),round(nrow(pt4)*0.8)); pt41  <- pt4[pt41s,]; pt42   <- pt4[-pt41s,] # T4
pt51s <- sample(1:nrow(pt5),round(nrow(pt5)*0.8)); pt51  <- pt5[pt51s,]; pt52   <- pt5[-pt51s,] # T5
pt61s <- sample(1:nrow(pt6),round(nrow(pt6)*0.8)); pt61  <- pt6[pt61s,]; pt62   <- pt6[-pt61s,] # T6



########################################
# Load features
Fmat <- matrix(NA,numPix,40) 
for(i in 1:numChunks){
  cc <- sprintf('%03d',i)
  cfile <- paste0(ckPheDir,'/chunk_feat_',cc,'.rda') 
  log <- try(load(cfile),silent=F)
  if (inherits(log, 'try-error')) next 
    
  if(i==numChunks){chunks <- c((chunk*(i-1)+1):numPix)
  }else{chunks <- c((chunk*(i-1)+1):(chunk*i))}
    
  chunkStart <- chunks[1];  chunkEnd <- chunks[length(chunks)]
    
  Fmat[chunkStart:chunkEnd,] <- f_mat[,(58*(yToDo-1)+c(1:40))]
  
 print(i)
}

Fmat11 <- Fmat[pt11[,2],]; Fmat12 <- Fmat[pt12[,2],]
Fmat21 <- Fmat[pt21[,2],]; Fmat22 <- Fmat[pt22[,2],]
Fmat31 <- Fmat[pt31[,2],]; Fmat32 <- Fmat[pt32[,2],]
Fmat41 <- Fmat[pt41[,2],]; Fmat42 <- Fmat[pt42[,2],]
Fmat51 <- Fmat[pt51[,2],]; Fmat52 <- Fmat[pt52[,2],]
Fmat61 <- Fmat[pt61[,2],]; Fmat62 <- Fmat[pt62[,2],]

lct <- c(rep(1,dim(Fmat11)[1]),rep(2,dim(Fmat21)[1]),rep(3,dim(Fmat31)[1]),
         rep(4,dim(Fmat41)[1]),rep(5,dim(Fmat51)[1]),rep(6,dim(Fmat61)[1]))



########################################
## Fit model
mat_input <- as.data.frame(rbind(Fmat11,Fmat21,Fmat31,Fmat41,Fmat51,Fmat61))

# # PCA to reduce dimensionallity
# pca <- princomp(~.,data=mat_input,cor=T,na.action=na.exclude)
# 
# plot(pca$sdev/sum(pca$sdev)*100)
# print(sum(as.numeric(pca$sdev/sum(pca$sdev)*100)[1:15])) # % of variance explained
# 
# # input
# input <- as.data.frame(cbind(lct,pca$scores[,1:15]))
# input <- na.omit(input)

# input
input <- as.data.frame(cbind(lct,mat_input))
input <- na.omit(input)


# Fit 
rf <- randomForest(input[,-1],as.factor(input[,1]),importance=T)
print(rf)


########################################
# pre_pc  <- predict(pca,as.data.frame(Fmat))
# cCase   <- complete.cases(pre_pc)
# pred_rf <- predict(rf,pre_pc[cCase,1:15])
# rf_val  <- matrix(NA,numPix,1)
# rf_val[cCase] <- pred_rf

cCase   <- complete.cases(Fmat)
pred_rf <- predict(rf,Fmat[cCase,])
rf_val  <- matrix(NA,numPix,1)
rf_val[cCase] <- pred_rf

pred_map <- setValues(imgBase,rf_val)


#######################################
# # Plot
# mycol <- c('#ff0004','#2535e1','#72dee0','#00ff22','#ff5cf7','#e1e300')
# plot(pred_map,col=mycol)

# Save
pheDir <- paste0(params$setup$outDir,strSite,'/output')
if (!dir.exists(pheDir)) {dir.create(pheDir)}

writeRaster(pred_map,filename=paste0(pheDir,'/map_re_1_',2017+yToDo,'.tif'),overwrite=TRUE)
save(rf,file=paste0(pheDir,'/rf_pc_',2017+yToDo,'.rda'))



# #######################################
# # Validation
# pheDir <- paste0(params$setup$outDir,strSite,'/output')
# pred_map <- rast(paste0(pheDir,'/map_re_',2017+yToDo,'.tif'))
# 
# pix_ref <- c(pt12[,2],pt22[,2],pt32[,2],pt42[,2],pt52[,2],pt62[,2])
# pts_ref <- values(pred_map)[pix_ref]
# pts_map <- c(rep(1,length(pt12[,2])),rep(2,length(pt22[,2])),rep(3,length(pt32[,2])),
#              rep(4,length(pt42[,2])),rep(5,length(pt52[,2])),rep(6,length(pt62[,2])))
# 
# print(sum(pts_ref[which(pts_map==1)]==pts_map[which(pts_map==1)])/length(pts_ref[which(pts_map==1)])*100)
# print(sum(pts_ref[which(pts_map==2)]==pts_map[which(pts_map==2)])/length(pts_ref[which(pts_map==2)])*100)
# print(sum(pts_ref[which(pts_map==3)]==pts_map[which(pts_map==3)])/length(pts_ref[which(pts_map==3)])*100)
# print(sum(pts_ref[which(pts_map==4)]==pts_map[which(pts_map==4)])/length(pts_ref[which(pts_map==4)])*100)
# print(sum(pts_ref[which(pts_map==5)]==pts_map[which(pts_map==5)])/length(pts_ref[which(pts_map==5)])*100)
# print(sum(pts_ref[which(pts_map==6)]==pts_map[which(pts_map==6)])/length(pts_ref[which(pts_map==6)])*100)


