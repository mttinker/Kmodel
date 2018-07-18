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
dfDU = read_xlsx("./Data/Depth_use_sum.xlsx")
dfHR = read_xlsx("./Data/Home Range Summaries.xlsx")
dfNgrids = read_xlsx("./Data/Count_grid_Strata.xlsx")
dfCarcs = read_xlsx("./Data/Shark_sum.xlsx")
# Import table of years when females began to reproduce in each hab section/Subpop
# (in years prior to this, a hab section was a male dominated area)
YearRep = read_xlsx("./Data/Repro_by_Subpop.xlsx"); YrFem = YearRep$YearRep

zeta=numeric(); rho = numeric()
zeta[1] = .01 # Intercept (nuiscence parameter)
zeta[2] = 1.3	# Shared hazards for all animals when density low, food abundant
zeta[3] = 0.5  # Additional hazards for males **ALLOW THIS TO BE FIT TO DATA?**
zeta[4] = 1.5  # Added density-dependent hazards (scaled by proportion of K)
rho[1] = 0.3 # Param for reproduction
rho[2] = 2  # Param for reproduction
rho[3] = -3 # Param for reproduction
rho[4] = 0.167 # Param for reproduction

dfGH$Count[is.na(dfGH$Count)]=0
dfGL$Count[is.na(dfGL$Count)]=0
dfGH$Depth = dfGH$Depth*-1 +1
dfGL$Depth = dfGL$Depth*-1 +1
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
LogInd = log(CountI+.5)
LogPup = log(CountP+.5)
logN0 = numeric()
for (i in 1:Npop){
  logN0[i] = (mean(LogInd[Psrv==i & Ysrv<7]))
}  

dfDU$Depth = dfDU$Depth*-1 +1
dfDU = dfDU[order(dfDU$Depth),]
Ndeps = max(dfDU$Depth)
NcellDep = dfDU$NumCells
CountDep = dfDU$SumOtts
RelDensDep = CountDep/(NcellDep*1.6)
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
Fmoves = c(rexp(100,1/rnorm(100,17.54,3.53)),
           rexp(100,1/rnorm(100,10.52,3.16)),
           rexp(100,1/rnorm(100,16.48,3.36)),
           rexp(300,1/rnorm(300,5.33,0.73)),
           rexp(300,1/rnorm(300,4.82,0.5)),
           rexp(300,1/rnorm(300,6.13,1.8)))  
Mmoves = c(rexp(100,1/rnorm(100,37.82,4.41)),
           rexp(100,1/rnorm(100,91.53,15.24)),
           rexp(100,1/rnorm(100,47.76,10.37)),
           rexp(300,1/rnorm(300,9.91,1.66)),
           rexp(300,1/rnorm(300,7.27,2.11)),
           rexp(300,1/rnorm(300,25.41,7.11)))  
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
# Summarize grid points for each hab section to use for fitting--------------------------------
Hsampsz = 24 # MUST BE EVEN
Lsampsz = 12 # MUST BE EVEN
# FIX THIS!!!
# NOTE : need to revise this!!! calculate values based on random sample of grid points?
GridcountH = matrix(data = 0,nrow = Hsampsz,ncol=Npop)
GridcountL = matrix(data = 0,nrow = Lsampsz,ncol=Npop)
#AreaGH = GridcountH
#PrpnAH = GridcountH
#AreaGL = GridcountL
#PrpnAL = GridcountL
DepH = GridcountH
KlpH = GridcountH
DshH = GridcountH
SubH = GridcountH
DepL = GridcountL
KlpL = GridcountL
DshL = GridcountL
SubL = GridcountL
SH1 = GridcountH
SH2 = GridcountH
SH3 = GridcountH
SH4 = GridcountH
SL1 = GridcountL
SL2 = GridcountL
SL3 = GridcountL
SL4 = GridcountL
NoEstu = GridcountH
AreaHD = numeric()
AreaLD = numeric()
NgrdH = numeric()
NgrdL = numeric()
#Kgrdtest = matrix(data = 0,nrow = Hsampsz+Lsampsz,ncol=Npop)
for (i in 1:Npop){
  if (i==5){  
    AreaHD[i] = .9999*Areas[i]
    AreaLD[i] = .0001*Areas[i]   
    NgrdH[i] = Hsampsz
    NgrdL[i] = 1    
    iiH = which(dfGH$Pop==i)
    iii = sample(iiH,Hsampsz)
    SubH[,i] = 4
    DepH[,i] = 3
    KlpH[,i] = 0
    DshH[,i] = 0   
    GridcountH[,i] = dfGH$Count[iii]
    # clust = ceiling(seq(1,Hsampsz,length.out = length(iiH)))
    # for (c in 1:max(clust)){
    #   AreaGH[c,i] = sum(clust==c)*.01
    #   PrpnAH[c,i] = sum(clust==c)/length(clust)
    #   GridcountH[c,i] = sum(ottH[clust==c])
    # }
    # AreaGH[1,i] = .01
    # PrpnAL[1,i] = 1
    SubL[,i] = 1
    DepL[,i] = 15
    KlpH[1,i] = 0
    DshH[1,i] = 0
  }else{
    NgrdH[i] = Hsampsz
    NgrdL[i] = Lsampsz
    AreaHD[i] = Areas[i]*(dfNgrids$HD[i]/(dfNgrids$HD[i]+dfNgrids$LD[i])) 
    AreaLD[i] = Areas[i]*(dfNgrids$LD[i]/(dfNgrids$HD[i]+dfNgrids$LD[i]))     
    iiH1 = which(dfGH$Pop==i & dfGH$Count>0)
    iiH0 = which(dfGH$Pop==i & dfGH$Count==0)
    prp = min(.9,5*length(iiH1)/(length(iiH1)+length(iiH0)))    
    selctPos = max(1,round(prp*Hsampsz))
    iii = c(sample(iiH1,selctPos),sample(iiH0,Hsampsz-selctPos))
    #iiH = which(dfGH$Pop==i) 
    #iii = sample(iiH,Hsampsz)
    GridcountH[,i] = dfGH$Count[iii]
    SubH[,i] = dfGH$SubstrCat[iii]
    DepH[,i] = dfGH$Depth[iii]
    KlpH[,i] = dfGH$Kelp[iii]
    DshH[,i] = dfGH$Distresid[iii]
    #
    iiL1 = which(dfGL$Pop==i & dfGL$Count>0)
    iiL0 = which(dfGL$Pop==i & dfGL$Count==0)
    prp = min(.9,10*length(iiL1)/(length(iiL1)+length(iiL0)))
    selctPos = max(1,round(prp*Lsampsz))
    iii = c(sample(iiL1,selctPos),sample(iiL0,Lsampsz-selctPos))    
    #iiL = which(dfGL$Pop==i)
    #iii = sample(iiL,Lsampsz)
    GridcountL[,i] = dfGL$Count[iii]
    SubL[,i] = dfGL$SubstrCat[iii]
    DepL[,i] = pmax(1,pmin(Ndeps,dfGL$Depth[iii]))
    KlpL[,i] = dfGL$Kelp[iii]
    DshL[,i] = dfGL$Distresid[iii]
    #     mat = as.matrix(cbind(dfGH$SubstrCat[iiH],dfGH$Depth[iiH],dfGH$Kelp[iiH],dfGH$Distresid[iiH]))
    # set.seed(20)
    # clustH = kmeans(mat,Hsampsz,nstart = 100)
    # clust = clustH$cluster
    # for (c in 1:max(clust)){
    #   AreaGH[c,i] = sum(clust==c)*.01
    #   PrpnAH[c,i] = sum(clust==c)/length(clust)
    #   SubH[c,i] = round(mean(mat[clust==c,1]))
    #   DepH[c,i] = round(mean(mat[clust==c,2]))
    #   KlpH[c,i] = mean(mat[clust==c,3])
    #   DshH[c,i] = mean(mat[clust==c,4])
    #   GridcountH[c,i] = sum(ottH[clust==c])
    # }
    # mat = as.matrix(cbind(dfGL$SubstrCat[iiL],dfGL$Depth[iiL],dfGL$Kelp[iiL],dfGL$Distresid[iiL]))
    # ottL = dfGL$Count[iiL]
    # set.seed(20)
    # clustL = kmeans(mat,Lsampsz,nstart = 100)
    # clust = clustL$cluster
    # for (c in 1:max(clust)){
    #   AreaGL[c,i] = sum(clust==c)*.01
    #   PrpnAL[c,i] = sum(clust==c)/length(clust)
    #   SubL[c,i] = round(mean(mat[clust==c,1]))
    #   DepL[c,i] = round(mean(mat[clust==c,2]))
    #   KlpL[c,i] = mean(mat[clust==c,3])
    #   DshL[c,i] = mean(mat[clust==c,4])
    #   GridcountL[c,i] = sum(ottL[clust==c])
    # }
  }
  NoEstu[which(SubH[,i]<4),i] = 1  
  SH1[which(SubH[,i]==1),i] = 1
  SH2[which(SubH[,i]==2),i] = 1
  SH3[which(SubH[,i]==3),i] = 1
  SH4[which(SubH[,i]==4),i] = 1
  SL1[which(SubL[,i]==1),i] = 1
  SL2[which(SubL[,i]==2),i] = 1
  SL3[which(SubL[,i]==3),i] = 1
  SL4[which(SubL[,i]==4),i] = 1
 # PrpnAH[,i] = PrpnAH[,i]/sum(PrpnAH[,i])
  # PrpnAL[,i] = PrpnAL[,i]/sum(PrpnAL[,i])
}
# Estimate proportion of cells in each depth bin for each pop
PrpnDH = matrix(nrow = Ndeps,ncol = Npop)
PrpnDL = matrix(nrow = Ndeps,ncol = Npop)
for (i in 1:Npop){
  for (d in c(1:61)){
    PrpnDH[d,i] = length(which(dfGH$Pop==i & dfGH$Depth==d))/length(which(dfGH$Pop==i))
    PrpnDL[d,i] = length(which(dfGL$Pop==i & dfGL$Depth==d))/length(which(dfGL$Pop==i & dfGL$Depth<=Ndeps))
  }
}
PrpnDL[,5] = 0; PrpnDL[15,5] = 1
#
# Calculate mean variable values for ALL gridcells in each strata -------------------
SH1all = numeric(length = length(dfGH$SubstrCat)); SH1all[which(dfGH$SubstrCat==1)]=1
SH2all = numeric(length = length(dfGH$SubstrCat)); SH2all[which(dfGH$SubstrCat==2)]=1
SH3all = numeric(length = length(dfGH$SubstrCat)); SH3all[which(dfGH$SubstrCat==3)]=1
SH4all = numeric(length = length(dfGH$SubstrCat)); SH4all[which(dfGH$SubstrCat==4)]=1
SL1all = numeric(length = length(dfGL$SubstrCat)); SL1all[which(dfGL$SubstrCat==1)]=1
SL2all = numeric(length = length(dfGL$SubstrCat)); SL2all[which(dfGL$SubstrCat==2)]=1
SL3all = numeric(length = length(dfGL$SubstrCat)); SL3all[which(dfGL$SubstrCat==3)]=1
SH1mn = numeric(); SH2mn = numeric(); SH3mn = numeric(); SH4mn = numeric()
SL1mn = numeric(); SL2mn = numeric(); SL3mn = numeric(); 
KlpHmn = numeric(); DshHmn = numeric(); KlpLmn = numeric(); DshLmn = numeric(); NoEstmn = numeric()
for (i in 1:Npop){
  SH1mn[i] = mean(SH1all[which(dfGH$Pop==i)])
  SH2mn[i] = mean(SH2all[which(dfGH$Pop==i)])
  SH3mn[i] = mean(SH3all[which(dfGH$Pop==i)])
  SH4mn[i] = mean(SH4all[which(dfGH$Pop==i)])
  KlpHmn[i] = mean(dfGH$Kelp[which(dfGH$Pop==i)])
  DshHmn[i] = mean(dfGH$Distresid[which(dfGH$Pop==i)])   
  NoEstmn[i]= mean(1-SH4all[which(dfGH$Pop==i)]) 
  if (i==5){
    SL1mn[i] = 1
    SL2mn[i] = 0
    SL3mn[i] = 0
    KlpLmn[i] = 0
    DshLmn[i] = 0
  }else{
    SL1mn[i] = mean(SL1all[which(dfGL$Pop==i)])
    SL2mn[i] = mean(SL2all[which(dfGL$Pop==i)])
    SL3mn[i] = mean(SL3all[which(dfGL$Pop==i)])  
    KlpLmn[i] = mean(dfGL$Kelp[which(dfGL$Pop==i)])
    DshLmn[i] = mean(dfGL$Distresid[which(dfGL$Pop==i)]) 
  }
}
# Use approx values of params to estimate "se" of log(Kdensity) values 
HG = Hsampsz; LG = Lsampsz
B1 = 8; B2 = 80; B3 = .9; B4 = -.5
B0 = numeric(); B0[1] = .7; B0[2] = .8; B0[3] = .9; B0[4] = 1.1;
DepEff <- pmax(.02,(1-((seq(1,61)-B1)/100)^2)^B2)
logDepEff = log(DepEff)
LogDepHD = numeric(); LogKmnH = numeric() ; 
LogDepLD = numeric(); LogKmnL = numeric() ; 
LogKavH = numeric() ; KavH = numeric() ; KmnH = numeric()
LogKavL = numeric() ; KavL = numeric() ; KmnL = numeric() 
seKH = numeric(); seKL = numeric();
for (i in 1:Npop){
  SH1ii = (SH1all[which(dfGH$Pop==i)])
  SH2ii = (SH2all[which(dfGH$Pop==i)])
  SH3ii = (SH3all[which(dfGH$Pop==i)])
  SH4ii = (SH4all[which(dfGH$Pop==i)])
  KlpHii = (dfGH$Kelp[which(dfGH$Pop==i)])
  DshHii = (dfGH$Distresid[which(dfGH$Pop==i)]) 
  DepHii = (dfGH$Depth[which(dfGH$Pop==i)])
  NoEstuii = (1-SH4all[which(dfGH$Pop==i)])
  if(i != 5){
    SL1ii = (SL1all[which(dfGL$Pop==i)])
    SL2ii = (SL2all[which(dfGL$Pop==i)])
    SL3ii = (SL3all[which(dfGL$Pop==i)])
    KlpLii = (dfGL$Kelp[which(dfGL$Pop==i)])
    DshLii = (dfGL$Distresid[which(dfGL$Pop==i)]) 
    DepLii = pmin(Ndeps,(dfGL$Depth[which(dfGL$Pop==i)]))  
  }
  # Calculate grid-cell expected log K and K vals
  LogKgHD = B0[1]*SH1ii+B0[2]*SH2ii+B0[3]*SH3ii+B0[4]*SH4ii+NoEstuii*(logDepEff[DepHii]+B3*KlpHii+B4*DshHii) 
  KgHD <- exp(LogKgHD)
  if(i != 5){
    LogKgLD = B0[1]*SL1ii+B0[2]*SL2ii+B0[3]*SL3ii+logDepEff[DepLii]+B3*KlpLii+B4*DshLii 
    KgLD <- exp(LogKgLD)    
  }else{
    LogKgLD = LogKgHD
    KgLD <- exp(LogKgLD)    
  }
  # Compute median-based std err (MAD) to avoid bias from large values
  if(mean(KgHD)>2){
    seKH[i] = 1.4826*median(abs(LogKgHD-median(LogKgHD)))
  }else{
    seKH[i] = sd(LogKgHD)
  }
  if (i == 5){
    seKL[i] = seKH[i] 
  }else{
    # seKL[i] = 1.4826*median(abs(LogKgLD-median(LogKgLD))) 
    seKL[i] = sd(LogKgLD)
  }
  LogDepHD[i] <- sum(PrpnDH[1:Ndeps,i]*logDepEff[1:Ndeps])
  LogKmnH[i] = B0[1]*SH1mn[i]+B0[2]*SH2mn[i]+B0[3]*SH3mn[i]+B0[4]*SH4mn[i]+NoEstmn[i]*(LogDepHD[i]+B3*KlpHmn[i]+B4*DshHmn[i])
  if(i != 5){
    LogDepLD[i] <- sum(PrpnDL[1:Ndeps,i]*logDepEff[1:Ndeps])
    LogKmnL[i] = B0[1]*SL1mn[i]+B0[2]*SL2mn[i]+B0[3]*SL3mn[i]+LogDepLD[i]+B3*KlpLmn[i]+B4*DshLmn[i]
  }else{
    LogDepLD[i] <- sum(PrpnDL[1:Ndeps,i]*logDepEff[1:Ndeps])
    LogKmnL[i] = LogKmnH[i]  
  }  
  print(" ")
  print(paste0("Subpop = ", i))
  print("Log mean K dens, HD, avg of grid vs mean pop")
  LogKavH[i] = mean(LogKgHD)
  print(LogKavH[i]) 
  print(LogKmnH[i])
  print("Est mean K dens, HD, avg of grid vs mean pop")
  KavH[i] = mean(KgHD)
  KmnH[i] = exp(LogKmnH[i]+(seKH[i]^2)/2)
  print(KavH[i]) 
  print(KmnH[i])
  print("Log mean K dens, LD, avg of grid vs mean pop")
  LogKavL[i] = mean(LogKgLD)
  print(LogKavL[i]) 
  print(LogKmnL[i])
  KavL[i] = mean(KgLD)
  KmnL[i] = exp(LogKmnL[i]+(seKL[i]^2)/2)
  print("Est mean K dens, LD, avg of grid vs mean pop")
  print(KavL[i]) 
  print(KmnL[i])
}
plot(KavH,KmnH)
lines(x = c(0,max(KavH)), y = c(0,max(KavH)))
plot(KavL,KmnL)
lines(x = c(0,max(KavL)), y = c(0,max(KavL)))
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
# dfGsmp = rbind(dfGH[iiH,],dfGL[iiL,]) 
# Ngrd = dim(dfGsmp)[1]
# hist(dfGsmp$SubstrCat)
save(dfHabSecs,dfC,dfGH,dfGL,disprobF,disprobM,distmat,
     DemComp,YearRep,Psrv,Ysrv,NCarcsF,NCarcsM,Pupratio,
     CountI,CountP,LogInd,LogPup,logN0,
     Areas,AreaHD,AreaLD,NgrdH,NgrdL,GridcountH,GridcountL,
     DepH,DepL,KlpH,KlpL,DshH,DshL,SubH,SubL,UB,Ncounts,
     CarcShrkF,PshkF,YshkF,CarcShrkM,PshkM,YshkM,
     TotCarcF,TotCarcM,dfDU,Ndeps,RelDensDep,HG,LG,
     PrpnDH,PrpnDL,NoEstu,NoEstmn,seKH,seKL,
     SH1,SH2,SH3,SH4,SL1,SL2,SL3,
     SH1mn,SH2mn,SH3mn,SH4mn,SL1mn,SL2mn,SL3mn,
     KlpHmn,DshHmn,KlpLmn,DshLmn,
     rho,zeta,Npop,Nyrs,Years,file = "Data_for_Kmod.rdata") 

