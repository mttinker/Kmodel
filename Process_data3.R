# Process census data, habitat data, and select random sample of grid cells
#
# **NOTE** FIX seKL and seKH so that it uses the clustered grid cells 
#
rm(list = ls())
library(readxl)
library(ggplot2)
library(gridExtra)
library(stats)
library(fitdistrplus)
library(extraDistr)
PlotTrends = 1
dfC = read_xlsx("./Data/Sum_Year_HabsecLim.xlsx")
dfD1 = read_xlsx("./Data/Sum_Pop_dep.xlsx")
dfD2 = read_xlsx("./Data/Sum_Pop_estuary.xlsx")
dfHabSecs = read_xlsx("./Data/HabSecsLim.xlsx")
# dfDU = read_xlsx("./Data/Depth_use_sum.xlsx")
# dfHR = read_xlsx("./Data/Home Range Summaries.xlsx")
# dfNgrids = read_xlsx("./Data/Count_grid_Strata.xlsx")
dfCarcs = read_xlsx("./Data/Shark_sum.xlsx")

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
LogInd = log(CountI+.001)
LogPup = log(CountP+.001)
logN0 = numeric()
for (i in 1:Npop){
  logN0[i] = (mean(LogInd[Psrv==i & Ysrv<7]))
}  
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
Depths = dfD1$Dep
Dists = dfD1$Dshore
dfDep = data.frame(Depth=Depths,Dist=Dists,logDist=log(Dists+1))
fitdep = nls(logDist ~ a*Depth^b+c, data = dfDep, start = list(a=25, b=.04, c = -20))
summary(fitdep)
deppars = coef(fitdep)      # model coefficients (means, slopes, intercepts)
depparCI = confint(fitdep)             # confidence intervals for parameters
pred = predict(fitdep, newdata=data.frame(Depth=1:60)) # predicted values for new observations
plot(dfDep$Depth,dfDep$logDist,xlab=c("Depth"),ylab=c("log(Distance from Shore)"))
lines(seq(1,60),pred, col="red")
pred = predict(fitdep)
resids = dfDep$logDist - pred
# hist(resids)
# plot(dfDep$Depth,resids)
a = deppars[1]; b = deppars[2]; c = deppars[3]
# Calc distshore residual valures for each point 
dfD1$Distresid = log(dfD1$Dshore+1) - (a*dfD1$Dep^b+c)
rm(a); rm(b); rm(c)
# 
# Create dispersal matrices for females and males
distmat = dist(0.5*dfHabSecs$Mean_ATOS); distmat = as.matrix(distmat)
# Correct Monterey Bay distances to allow cross-bay movements?
# distmat[3,7] = 35; distmat[7,3] = 35; 
# distmat[4,7] = 30; distmat[7,4] = 30; 
# distmat[5,7] = 25; distmat[7,5] = 25; 
# distmat[6,7] = 15; distmat[7,6] = 15; 
# Fmoves = numeric()
# Mmoves = numeric()
# for (i in 1:dim(dfHR)[1]){
#   if(dfHR$Sex[i]=="F"){
#     Fmoves = c(Fmoves,rbeta(100,1,4)*dfHR$Range_span_km[i])
#   }else{
#     Mmoves = c(Mmoves,rbeta(100,1,4)*dfHR$Range_span_km[i])
#   }
# }
Fmoves = c(rexp(100,1/pmax(.001,rnorm(100,17.54,3.53))),
           rexp(100,1/pmax(.001,rnorm(100,10.52,3.16))),
           rexp(100,1/pmax(.001,rnorm(100,16.48,3.36))),
           rexp(500,1/pmax(.001,rnorm(500,5.33,0.73))),
           rexp(500,1/pmax(.001,rnorm(500,4.82,0.5))),
           rexp(500,1/pmax(.001,rnorm(500,6.13,1.8))))  
Mmoves = c(rexp(100,1/pmax(.001,rnorm(100,37.82,4.41))),
           rexp(100,1/pmax(.001,rnorm(100,91.53,15.24))),
           rexp(100,1/pmax(.001,rnorm(100,47.76,10.37))),
           rexp(500,1/pmax(.001,rnorm(500,9.91,1.66))),
           rexp(500,1/pmax(.001,rnorm(500,7.27,2.11))),
           rexp(500,1/pmax(.001,rnorm(500,25.41,7.11))))  
dispfitF = fitdist(Fmoves,distr = 'weibull')
summary(dispfitF); # plot(dispfitF)
dispfitM = fitdist(Mmoves,distr = 'weibull')
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
# Process shark data-----------------------------------------------
CarcShrkF = numeric()
CarcShrkM = numeric()
TotCarcF = numeric()
TotCarcM = numeric()
PshkF = numeric()
YshkF = numeric()
PshkM = numeric()
YshkM = numeric()
NCarcsF = 0
NCarcsM = 0
for (i in 1:Nyrs){
  for (j in 1:Npop){
    ii = which(dfCarcs$SEX=="F" & dfCarcs$Pop==j & dfCarcs$Year==Years[i])
    if(length(ii)>0){
      NCarcsF = NCarcsF+1
      PshkF = c(PshkF,j)
      YshkF = c(YshkF,i)
      TotCarcF = c(TotCarcF,length(ii))
      CarcShrkF = c(CarcShrkF,sum(dfCarcs$Shark[ii]))
    }
    ii = which(dfCarcs$SEX=="M" & dfCarcs$Pop==j & dfCarcs$Year==Years[i])
    if(length(ii)>0){
      NCarcsM = NCarcsM+1
      PshkM = c(PshkM,j)
      YshkM = c(YshkM,i)
      TotCarcM = c(TotCarcM,length(ii))
      CarcShrkM = c(CarcShrkM,sum(dfCarcs$Shark[ii]))
    }    
  }
}
# Process habitat variables for each Pop: make depth bin matrices------
NbinAvg = 27
Nbin = numeric(length=Npop)
AreaB = matrix(nrow = NbinAvg, ncol = Npop)
Dep = matrix(nrow = NbinAvg, ncol = Npop)
Phrd = matrix(nrow = NbinAvg, ncol = Npop)
Pklp = matrix(nrow = NbinAvg, ncol = Npop)
Dsh = matrix(nrow = NbinAvg, ncol = Npop)
Sumcount = matrix(nrow = NbinAvg, ncol = Npop)
for (i in 1:Npop){
  ii = which(dfD2$Pop==i)
  if(length(ii)>0){
    AreaB[1,i] = dfD2$Ncells[ii]*0.01
    Sumcount[1,i] = dfD2$SumInd[ii]
  }else{
    AreaB[1,i] = 0
    Sumcount[1,i] = 0
  }
  ii = which(dfD1$Pop==i)
  Nbin[i] = length(ii)+1
  AreaB[2:Nbin[i],i] = dfD1$Ncells[ii]*.01
  Dep[2:Nbin[i],i] = dfD1$Dep[ii]
  Phrd[2:Nbin[i],i] = dfD1$PpnHard[ii]
  Pklp[2:Nbin[i],i] = dfD1$Kelp[ii]
  Dsh[2:Nbin[i],i] = dfD1$Distresid[ii]
  Sumcount[2:Nbin[i],i] = dfD1$SumInd[ii]
}
AreaHD = numeric()
for (i in 1:Npop){
  sumArea = sum(AreaB[,i],na.rm = TRUE)
  AreaB[,i] = Areas[i]*(AreaB[,i]/sumArea)
  AreaHD[i] = sum(AreaB[1:23,i],na.rm = TRUE)
}  

# Estimate starting point for K values-------------------------
Kguess = matrix(nrow=6,ncol=Npop)
seKg = numeric()
for (i in 1:Npop){
  if(i==5){
    iii = which(dfC$Pop==i & dfC$Year>=2012 & dfC$Year<=2017)
  }else{
    iii = which(dfC$Pop==i & dfC$Year>=2005 & dfC$Year<=2010)
  }
  Kguess[1:6,i] = log(dfC$Otts[iii]/AreaHD[i])
  seKg[i] = sd(Kguess[1:6,i])
}
seKg = 1.5*as.numeric(quantile(seKg,.75))
NKguess=dim(Kguess)[1]
# dfGsmp = rbind(dfGH[iiH,],dfGL[iiL,]) 
# Ngrd = dim(dfGsmp)[1]
# hist(dfGsmp$SubstrCat)
save(dfHabSecs,dfC,disprobF,disprobM,distmat,
     DemComp,Psrv,Ysrv,NCarcsF,NCarcsM,Pupratio,
     CountI,CountP,LogInd,LogPup,logN0,Areas,
     Nbin,AreaB,Dep,Phrd,Pklp,Dsh,UB,Ncounts,
     CarcShrkF,PshkF,YshkF,CarcShrkM,PshkM,YshkM,
     TotCarcF,TotCarcM,Sumcount,Kguess,seKg,NKguess,
     rho,zeta,Npop,Nyrs,Years,AreaHD,
     file = "Data_for_Kmod3.rdata") 

