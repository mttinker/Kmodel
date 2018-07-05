# Process census data, habitat data, and select random sample of grid cells
library(readxl)
library(ggplot2)
library(gridExtra)
library(stats)
library(fitdistrplus)
library(extraDistr)
PlotTrends = 1
dfC = read_xlsx("./Data/Sum_Year_HabsecLim.xlsx")
dfGH = read_xlsx("./Data/Sum_by_Grid_HD.xlsx")
dfGL = read_xlsx("./Data/Sum_by_Grid_LD.xlsx")
dfHabSecs = read_xlsx("./Data/HabSecsLim.xlsx")
dfdsp = read_xlsx("./Data/ottmovesdat.xlsx")
dfNgrids = read_xlsx("./Data/Count_grid_Strata.xlsx")
dfCarcs = read_xlsx("./Data/Shark_sum.xlsx")
# Import table of years when females began to reproduce in each hab section/Subpop
# (in years prior to this, a hab section was a male dominated area)
YearRep = read_xlsx("./Data/Repro_by_Subpop.xlsx"); YrFem = YearRep$YearRep
dfGH$Count[is.na(dfGH$Count)]=0
dfGL$Count[is.na(dfGL$Count)]=0
dfGH$Depth = dfGH$Depth*-1
dfGL$Depth = dfGL$Depth*-1
dfGH$Pop = dfGH$Hab_Sec-2
dfGL$Pop = dfGL$Hab_Sec-2
dfC$Pop = dfC$HabSecLim-2
Ncounts = dim(dfC)[1]
Npop = max(dfC$Pop)
Nyrs = max(dfC$Year)-min(dfC$Year)+1
Years = seq(min(dfC$Year),max(dfC$Year))
dfC$yearN = dfC$Year-min(dfC$Year)+1
Areas = dfHabSecs$AREA_KM2
Psrv = dfC$Pop
Ysrv = dfC$yearN
CountI = dfC$Otts
CountP = dfC$Pups
#
Pupratio = matrix(data = 0, nrow=Npop,ncol=Nyrs)
if (PlotTrends==1){
  par(mfrow=c(5,4))
  for (i in 1:Npop){
    ii = which(dfC$Pop==i)
    plot(dfC$Year[ii],dfC$PupRatio[ii],xlab=paste0("Year, SubPop ",i),
         ylab="PupRatio")
    fit1<-smooth.spline(dfC$Year[ii],dfC$PupRatio[ii],df=7) #16 degrees of freedom
    Pupratio[i,] = pmax(.01,predict(fit1,1983:2017)$y)
    lines(fit1,col="red",lwd=2)
  }
  # Import table of years when females began to reproduce in each hab section/Subpop
  # (in years prior to this, a hab section was a male dominated area)
  # Plot trends
  par(mfrow=c(5,4))
  for (i in 1:Npop){
    ii = which(dfC$Pop==i)
    plot(dfC$Year[ii],dfC$Density[ii],xlab=paste0("Year, SubPop ",i),
         ylab="Density")
    fit1<-smooth.spline(dfC$Year[ii],dfC$Density[ii],df=10) #16 degrees of freedom
    lines(fit1,col="red",lwd=2)
  }
  par(mfrow=c(1,1))
}
# Distance from Shore Function (creates DishShore residuals)
# (note: modal depth ~ 10)
sampsz = 2000
iix = sample(dim(dfGH)[1],sampsz)
iiy = sample(dim(dfGL)[1],sampsz)
Depths = c(as.numeric(dfGH$Depth[iix]),as.numeric(dfGL$Depth[iiy])) 
Dists =  c(as.numeric(dfGH$DistShore[iix]),as.numeric(dfGL$DistShore[iiy])) 
dfDep = data.frame(Depth=Depths,Dist=Dists,logDist=log(Dists+1))
fitdep = nls(logDist ~ a*Depth^b+c, data = dfDep, start = list(a=1, b=1, c = 2))
summary(fitdep)
deppars = coef(fitdep)      # model coefficients (means, slopes, intercepts)
depparCI = confint(fitdep)             # confidence intervals for parameters
pred = predict(fitdep, newdata=data.frame(Depth=0:70)) # predicted values for new observations
iii = sample(2*sampsz,1000)
# plot(dfDep$Depth[iii],dfDep$logDist[iii],xlab=c("Depth"),ylab=c("log(Distance from Shore)"))
# lines(seq(0,70),pred, col="red")
pred = predict(fitdep)
resids = dfDep$logDist - pred
# hist(resids)
# plot(dfDep$Depth,resids)
a = deppars[1]; b = deppars[2]; c = deppars[3]
# Calc distshore residual valures for each point 
dfGH$Distresid = log(dfGH$DistShore+1) - (a*dfGH$Depth^b+c)
dfGL$Distresid = log(dfGL$DistShore+1) - (a*dfGL$Depth^b+c)
rm(a); rm(b); rm(c)
# 
# Create dispersal matrices for females and males
distmat = dist(0.5*dfHabSecs$Mean_ATOS); distmat = as.matrix(distmat)
# Correct Monterey Bay distances to allow cross-bay movements
distmat[3,7] = 35; distmat[7,3] = 35; 
distmat[4,7] = 30; distmat[7,4] = 30; 
distmat[5,7] = 25; distmat[7,5] = 25; 
distmat[6,7] = 15; distmat[7,6] = 15; 
dispfitF = fitdist(dfdsp$move[dfdsp$Sex==0],distr = 'weibull')
summary(dispfitF); # plot(dispfitF)
dispfitM = fitdist(dfdsp$move[dfdsp$Sex==1],distr = 'weibull')
summary(dispfitM); # plot(dispfitM)
disprobF = matrix(data=0, nrow = Npop, ncol = Npop)
disprobM = matrix(data=0, nrow = Npop, ncol = Npop)
for (i in 1:Npop){
  disprobF[,i] = dweibull(distmat[,i],dispfitF$estimate[1],dispfitF$estimate[2])
  disprobF[i,i] = 0
  disprobF[,i] = disprobF[,i]/sum(disprobF[,i])
  probstay = pweibull(min(distmat[-i,i])/2,dispfitF$estimate[1],dispfitF$estimate[2])
  disprobF[,i] = disprobF[,i]*(1-probstay)
  disprobF[i,i] = probstay
  disprobM[,i] = dweibull(distmat[,i],dispfitM$estimate[1],dispfitM$estimate[2])
  disprobM[i,i] = 0
  disprobM[,i] = disprobM[,i]/sum(disprobM[,i])
  probstay = pweibull(min(distmat[-i,i])/2,dispfitM$estimate[1],dispfitM$estimate[2])
  disprobM[,i] = disprobM[,i]*(1-probstay)
  disprobM[i,i] = probstay  
}
# Create matrices to keep track of male areas (reduced female movements to male areas)
# and "urchin boom" impacts (2013-17, central coast)
DemComp = matrix(data = 0, nrow=Npop,ncol=Nyrs)
UB = matrix(data = 0, nrow=Npop,ncol=Nyrs)
for (i in 1:Npop){
  for (t in 1:Nyrs){
    if(Years[t]>=2013 & (dfHabSecs$Pop[i]>=7 &  dfHabSecs$Pop[i]<=9)){
      UB[i,t] = 1
    }       
    if(t==1){
      if(Pupratio[i,min(Nyrs,t+3)]<0.05){
        DemComp[i,t]=1
      }else if(Pupratio[i,min(Nyrs,t+3)]>0.1){
        DemComp[i,t]=3
      }else{
        DemComp[i,t]=2
      }       
    }else{
      if(Pupratio[i,min(Nyrs,t+3)]<0.05 & DemComp[i,(t-1)]==1){
        DemComp[i,t]=1
      }else if(Pupratio[i,min(Nyrs,t+3)]>0.1){
        DemComp[i,t]=3
      }else{
        DemComp[i,t]=2
      } 
    }
  }
}
#
# Sample random grid points for each hab section to use for fitting
Hsampsz = 300
Lsampsz = 100
GridcountH = matrix(data = 0,nrow = Hsampsz,ncol=Npop)
GridcountL = matrix(data = 0,nrow = Lsampsz,ncol=Npop)
DepH = GridcountH
KlpH = GridcountH
DshH = GridcountH
SubH = GridcountH
DepL = GridcountL
KlpL = GridcountL
DshL = GridcountL
SubL = GridcountL
AreaHD = numeric()
AreaLD = numeric()
NgrdH = numeric()
NgrdL = numeric()
for (i in 1:Npop){
  if (i==5){
    iiH = sample(which(dfGH$Pop==i),200) 
    AreaHD[i] = Areas[i]
    AreaLD[i] = 0
    NgrdH[i] = length(iiH)
    NgrdL[i] = 1
    SubL[1,i] = 1
    DepL[1,i] = 5
    GridcountH[1:NgrdH[i],i] = dfGH$Count[iiH]
    DepH[1:NgrdH[i],i] = dfGH$Depth[iiH]
    KlpH[1:NgrdH[i],i] = dfGH$Kelp[iiH]
    DshH[1:NgrdH[i],i] = dfGH$Distresid[iiH]
    SubH[1:NgrdH[i],i] = dfGH$SubstrCat[iiH]
  }else{
    NgrdH[i] = Hsampsz
    NgrdL[i] = Lsampsz
    iiH = sample(which(dfGH$Pop==i),Hsampsz) 
    iiL = sample(which(dfGL$Pop==i & dfGL$Depth<=60),Lsampsz)
    AreaHD[i] = Areas[i]*(dfNgrids$HD[i]/(dfNgrids$HD[i]+dfNgrids$LD[i])) 
    AreaLD[i] = Areas[i]*(dfNgrids$LD[i]/(dfNgrids$HD[i]+dfNgrids$LD[i])) 
    GridcountH[1:NgrdH[i],i] = dfGH$Count[iiH]
    GridcountL[1:NgrdL[i],i] = dfGL$Count[iiL]
    DepH[1:NgrdH[i],i] = dfGH$Depth[iiH]
    DepL[1:NgrdL[i],i] = dfGL$Depth[iiL]
    KlpH[1:NgrdH[i],i] = dfGH$Kelp[iiH]
    KlpL[1:NgrdL[i],i] = dfGL$Kelp[iiL]
    DshH[1:NgrdH[i],i] = dfGH$Distresid[iiH]
    DshL[1:NgrdL[i],i] = dfGL$Distresid[iiL]
    SubH[1:NgrdH[i],i] = dfGH$SubstrCat[iiH]
    SubL[1:NgrdL[i],i] = dfGL$SubstrCat[iiL]
  }
}
zeta=numeric(); rho = numeric()
zeta[1] = .01 # Intercept (nuiscence parameter)
zeta[2] = 1.3	# Shared hazards for all animals when density low, food abundant
zeta[3] = 0.5  # Additional hazards for males **ALLOW THIS TO BE FIT TO DATA?**
zeta[4] = 1.5  # Added density-dependent hazards (scaled by proportion of K)
rho[1] = 0.3 # Param for reproduction
rho[2] = 2  # Param for reproduction
rho[3] = -3 # Param for reproduction
rho[4] = 0.167 # Param for reproduction
# 
# Process shark data
ii = which(dfCarcs$SEX=="F")
NCarcsF = length(ii)
CarcShrkF = dfCarcs$Shark[ii]
PshkF = dfCarcs$Pop[ii]
YshkF = dfCarcs$Year[ii] - min(Years) + 1
ii = which(dfCarcs$SEX=="M")
NCarcsM = length(ii)
CarcShrkM = dfCarcs$Shark[ii]
PshkM = dfCarcs$Pop[ii]
YshkM = dfCarcs$Year[ii] - min(Years) + 1
# dfGsmp = rbind(dfGH[iiH,],dfGL[iiL,]) 
# Ngrd = dim(dfGsmp)[1]
# hist(dfGsmp$SubstrCat)
save(dfHabSecs,dfC,dfGH,dfGL,disprobF,disprobM,distmat,DemComp,
     YearRep,Psrv,Ysrv,CountI,CountP,NCarcsF,NCarcsM,Pupratio,
     Areas,AreaHD,AreaLD,NgrdH,NgrdL,GridcountH,GridcountL,
     DepH,DepL,KlpH,KlpL,DshH,DshL,SubH,SubL,UB,Ncounts,
     CarcShrkF,PshkF,YshkF,CarcShrkM,PshkM,YshkM,
     rho,zeta,Npop,Nyrs,Years,file = "Data_for_Kmod.rdata") 
