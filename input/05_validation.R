library(sp)
library(raster)
library(terra)
library(sf)
library(rjson)



#######################################
# Area estimation
map_pmf  <- rast('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/map_re_2020.tif')

map_gmw30 <- vect('/projectnb/modislc/users/mkmoon/mangrove/others/gmw/GMW_v3_2020/gmw_v3_2020_vec.shp')  
map_hgmf  <- vect('/projectnb/modislc/users/mkmoon/mangrove/others/GlobalMangrove2020/GlobalMangrove2020.shp')  

map_gmw30 <- project(map_gmw30,map_pmf)
map_hgmf  <- project(map_hgmf,map_pmf)

map_gmw30 <- crop(map_gmw30,map_pmf)
map_hgmf  <- crop(map_hgmf,map_pmf)


plot(map_pmf)
plot(map_gmw30,add=T,border='red')
plot(map_hgmf,add=T,border='blue')

plot(map_gmw30)
plot(map_hgmf)

# rasterize
rast_1 <- rasterize(map_gmw30,map_pmf)
plot(rast_1)
rast_2 <- rasterize(map_hgmf,map_pmf)
plot(rast_2)

print(sum(values(map_pmf)==1,na.rm=T)*9)
print(sum(values(rast_1)==1,na.rm=T)*9)
print(sum(values(rast_2)==1,na.rm=T)*9)


# changes 
map_1  <- rast('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/map_re_2018.tif')
map_2  <- rast('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/map_re_2019.tif')
map_3  <- rast('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/map_re_2020.tif')
map_4  <- rast('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/map_re_2021.tif')
map_5  <- rast('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/map_re_2022.tif')
map_6  <- rast('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/map_re_2023.tif')

dat <- matrix(c(sum(values(map_1)==1,na.rm=T)*9,
                sum(values(map_2)==1,na.rm=T)*9,
                sum(values(map_3)==1,na.rm=T)*9,
                sum(values(map_4)==1,na.rm=T)*9,
                sum(values(map_5)==1,na.rm=T)*9,
                sum(values(map_6)==1,na.rm=T)*9),6,1)
plot(dat)
summary(dat)
sd(dat)


#######################################
load('/projectnb/modislc/users/mkmoon/mangrove/Img_cliped/sanur/output/rf_2020.rda')

vari <- c('Blue Min.','Blue Max.','Blue Medi.','Blue Area','Blue Amp.',
          'Green Min.','Green Max.','Green Medi.','Green Area','Green Amp.',
          'Red Min.','Red Max.','Red Medi.','Red Area','Red Amp.',
          'NIR Min.','NIR Max.','NIR Medi.','NIR Area','NIR Amp.',
          'NDVI Min.','NDVI Max.','NDVI Medi.','NDVI Area','NDVI Amp.',
          'GNDVI Min.','GNDVI Max.','GNDVI Medi.','GNDVI Area','GNDVI Amp.',
          'NIRv Min.','NIRv Max.','NIRv Medi.','NIRv Area','NIRv Amp.',
          'EVI2 Min.','EVI2 Max.','EVI2 Medi.','EVI2 Area','EVI2 Amp.')

aa <- rf$importance
od <- order(aa[,8],decreasing=F)


library(RColorBrewer)
mycol <- brewer.pal(8,'Set1')
mycol <- c(mycol[5],mycol[5],mycol[5],mycol[5],mycol[5],
           mycol[6],mycol[6],mycol[6],mycol[6],mycol[6],
           mycol[7],mycol[7],mycol[7],mycol[7],mycol[7],
           mycol[8],mycol[8],mycol[8],mycol[8],mycol[8],
           mycol[2],mycol[2],mycol[2],mycol[2],mycol[2],
           mycol[3],mycol[3],mycol[3],mycol[3],mycol[3],
           mycol[4],mycol[4],mycol[4],mycol[4],mycol[4],
           mycol[1],mycol[1],mycol[1],mycol[1],mycol[1])

setwd('/projectnb/modislc/users/mkmoon/mangrove/figures/')
png(filename='importance.png',width=8.5,height=9,unit='in',res=300)

par(oma=c(1,4.5,1,1),mar=c(4,4,1,1),mgp=c(2.5,1,0))
barplot(aa[od[26:40],8],names.arg=vari[od[26:40]],horiz=T,col=mycol[od[26:40]],las=1,cex.axis=1.8,cex.names=1.5,
        xlab='Mean Decrease Gini',cex.lab=1.8)

dev.off()
