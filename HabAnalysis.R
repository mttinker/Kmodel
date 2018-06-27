library(rgdal)
library(raster)
library(readxl)
# Determine modal depth
dfDep = read_xlsx("./Data/Counts_by_Dep.xlsx")
dfDep$Dep = -1*dfDep$Depth
dfDep = dfDep[dfDep$Dep<=60,]
dfDep$Dep2 = dfDep$Dep^2
dfDep$Dep3 = dfDep$Dep^3
dfDep$Count = 4*dfDep$Freq
dfDep$Failures = dfDep$Count-dfDep$SumCounts
fit1 = lm(RelDens~poly(Dep,5),weights = Freq, data = dfDep)
summary(fit1)
pred <- predict(fit1,newdat = data.frame(Dep=seq(0,35)))
plot(dfDep$Dep, dfDep$RelDens, type="o", xlab="Depth", ylab="Freq_Occurence",xlim=c(0,35))
lines(seq(0,35), pred, lwd=2, col="red")
iii = which(pred==max(pred))
ModalDep = iii
#
# Use this to smooth and classify habitat variables from grid cell centroid layer
dfgrd = read.csv("./Data/Ottrange_Grid_CC.csv") # Load csv pointfile of grid centroids
# Add coordinates
dfgrd$X = (dfgrd$left+dfgrd$right)/2
dfgrd$Y = (dfgrd$top+dfgrd$bottom)/2
# Create California Teale Albers Projection object
TA <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +datum=NAD83 +units=m +ellps=GRS80 +towgs84=0,0,0")
# Create base raster grid for Central California coast
xmn=min(dfgrd$X); xmx=max(dfgrd$X); ymn=min(dfgrd$Y); ymx=max(dfgrd$Y)
baserast = raster(xmn=xmn, xmx=xmx, ymn=ymn, ymx=ymx, 
                   nrows=(xmx-xmn)/100, ncols=(ymx-ymn)/100, crs=TA)
# Ratserize dfgrid Substrate layer (-1 = estuary, 0 = soft, 1 = hard)
grdSub = rasterize(dfgrd[, c('X', 'Y')], baserast, dfgrd[, 'SubstrN'], fun=mean)
extMont = extent(c(-180000, -165000, -165000, -148000))
extElk = extent(c(-165000, -150000, -135000, -125000))
pal <- colorRampPalette(c("orange","purple"))
plot(grdSub,ext = extMont,col = pal(10))
plot(grdSub,ext = extElk,col = pal(10))

# re-sample using 500m moving window average (circle), to make smoothed substrate layer
fw <- focalWeight(grdSub, 250, type='circle') 
grdSubSm<-focal(grdSub, w=fw,fun=mean , na.rm=TRUE) 
plot(grdSubSm,ext = extMont,col = pal(10))
plot(grdSubSm,ext = extElk,col = pal(10))
sp <- SpatialPoints(data.frame(X = dfgrd$X,Y=dfgrd$Y))
Subsm = raster::extract(grdSubSm,sp,method='bilinear')
dfgrd$SubstrNsm = pmax(-1,pmin(1,Subsm*20))
ii = which(is.na(Subsm))
dfgrd$SubstrNsm[ii] = dfgrd$SubstrN[ii]
#
# evaluate results and classify as estuary, soft, mixed, hard
hist(dfgrd$SubstrNsm[dfgrd$SubstrN==-1],main="Estuary")
hist(dfgrd$SubstrNsm[dfgrd$SubstrN==0],main="Soft")
hist(dfgrd$SubstrNsm[dfgrd$SubstrN==1],main="Hard") # &dfgrd$SubstrNsm<.9
hist(dfgrd$SubstrNsm[dfgrd$SubstrNsm<.9],main="All")
dfgrd$SubstrCat = dfgrd$SubstrN
ii = which(dfgrd$SubstrNsm < -.5); dfgrd$SubstrCat[ii] = 4 # Estuary
ii = which(dfgrd$SubstrNsm >= .75); dfgrd$SubstrCat[ii] = 3 # Hard
ii = which(dfgrd$SubstrNsm >= -.5 & dfgrd$SubstrNsm <=.1); dfgrd$SubstrCat[ii] = 1 # Soft
ii = which(dfgrd$SubstrNsm > .1 & dfgrd$SubstrNsm <.75); dfgrd$SubstrCat[ii] = 2 # Mixed
# Map results of classification:
grdSubcl = rasterize(dfgrd[, c('X', 'Y')], baserast, dfgrd[, 'SubstrCat'], fun=mean)
plot(grdSubcl,ext = extMont,col = pal(3))
plot(grdSubcl,ext = extElk,col = pal(3))
write.csv(dfgrd,file="./Data/Ottrange_Grid_CC_class.csv",row.names = FALSE)



