# Script to estimate proportion K and Est N by Area
# Load results file -----------------------------------------------
load("../Results/FitK_Results_mod3.rdata")
require(runjags)
require(parallel)
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(dplyr)
require(reshape2)
require(stats)
require(xlsx)
Nsimsamp = sample(Nsims,1000,replace = TRUE)

# Stats by sub-region/Area -------------------------------------------
# Year and Area defs
Yr1 = 16
Yr2 = 30
YrsN = seq(Yr1,Yr2+1)
Yrs = seq(1998,2012)

AreaNames = c('NC','MB','RC','SC','PC','SN')
NpopAreas = c(3,3,5,3,5)
Areadefs = matrix(c(1, 1, 3,
                       2, 4, 6,
                       3, 7, 11,
                       4, 12,14,
                       5, 15,19),
            nrow = 5, ncol = 3, byrow = TRUE)
# Sum N estimates by year
EstN = matrix(nrow = 5,ncol = 16)
EstK = numeric(length = 5)
EstpK = matrix(nrow = 6,ncol = 15)
Estlam = matrix(nrow = 6,ncol = 15)
for (i in 1:5){
  t = 1
  tt = YrsN[t]
  ppi = seq(Areadefs[i,2],Areadefs[i,3])
  nn = numeric()
  kk = numeric()
  for (pp in ppi){
    nn = sum(nn,sumstats[which(startsWith(vn,paste0("N[",tt,",")) 
                               & endsWith(vn, paste0(",",pp,"]"))),4])
    kk = sum(kk,sumstats[which(vn==paste0("K[",pp,"]")),2])
  }
  EstN[i,1] = nn
  EstK[i] = kk
  for (t in 2:16) {
    tt = YrsN[t]
    nn = numeric()
    kk = numeric()
    for (pp in ppi){
      nn = sum(nn,sumstats[which(startsWith(vn,paste0("N[",tt,",")) 
                                 & endsWith(vn, paste0(",",pp,"]"))),4])
    }
    EstN[i,t] = nn
    Estlam[i,t-1] = EstN[i,t]/EstN[i,t-1]
    EstpK[i,t-1] = EstN[i,t-1]/EstK[i]
  }
}
# Load San Nicolas data
EstpK[6,] = c(0.071705718,	0.08249284,	0.092324095,	0.102802347,		
  0.112221213,  0.126051316,	0.134007995,	0.129348883,	0.132375508,	
  0.125599703, 0.136391973,	0.143389711,	0.165229953,	0.189074129,
  0.215803276)

Estlam[6,] = c(1.157955597,	1.127338868,	1.119178664,	1.09970169,	
               1.12960738,	1.065785248,	0.972006154,	1.030732016,
               0.953574425,	1.091887339,	1.055486777,	1.160408094,	
               1.151090735,	1.147817728,	1.088185806)
# Smooth
tmp1 = EstpK; tmp2 = Estlam
for (i in 1:6){
  tmp1[i,] = predict(loess(EstpK[i,]~seq(1,15), span = 0.7))
  tmp2[i,] = predict(loess(Estlam[i,]~seq(1,15), span = 0.7))
  plot(Yrs,EstpK[i,], main = paste0('Est pK, ', AreaNames[i])) 
  lines(Yrs,tmp1[i,])
  plot(Yrs,Estlam[i,], main = paste0('Est lam, ', AreaNames[i])) 
  lines(Yrs,tmp2[i,])
}
EstpK = tmp1
Estlam = tmp2

df_pK = data.frame(Area = AreaNames, AreaN = seq(1,6), Y1998=EstpK[,1],Y1999=EstpK[,2],
                   Y2000=EstpK[,3],Y2001=EstpK[,4],Y2002=EstpK[,5],
                   Y2003=EstpK[,6],Y2004=EstpK[,7],Y2005=EstpK[,8],
                   Y2006=EstpK[,9],Y2007=EstpK[,10],Y2008=EstpK[,11],
                   Y2009=EstpK[,12],Y2010=EstpK[,13],Y2011=EstpK[,14],
                   Y2012=EstpK[,15])
write.xlsx(df_pK,'../Results/Est_pK_Areas.xlsx',row.names = FALSE)

df_Lam = data.frame(Area = AreaNames, AreaN = seq(1,6), Y1998=Estlam[,1],Y1999=Estlam[,2],
                   Y2000=Estlam[,3],Y2001=Estlam[,4],Y2002=Estlam[,5],
                   Y2003=Estlam[,6],Y2004=Estlam[,7],Y2005=Estlam[,8],
                   Y2006=Estlam[,9],Y2007=Estlam[,10],Y2008=Estlam[,11],
                   Y2009=Estlam[,12],Y2010=Estlam[,13],Y2011=Estlam[,14],
                   Y2012=Estlam[,15])
write.xlsx(df_Lam,'../Results/Est_Lam_Areas.xlsx',row.names = FALSE)

