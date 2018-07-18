# Script to fit a state-space population model with density-dependence
# to various data sets from California, including annual census counts,
# average abundance at regularly spaced 100m grid cells along the coast,
# And stranding data on cause of death ("Shark" vs "other").
# Load libraries -------------------------------------------------------------
#require(readxl)
#require(stats)
#require(MASS)
#require(fitdistrplus)
require(runjags)
require(parallel)
require(gtools)
#require(lme4)
require(lattice)
require(coda)
# require(boot)
# require(rjags)
# require(doParallel)
# require(ggplot2)
# Load data ------------------------------------------------------------------
# source("Process_data.R") # NOTE: need to run Process_data script first
rm(list = c())
load("Data_for_Kmod.rdata")
# Starting proportion females for male areas, transition and female areas:
propfem = c(0.05,0.1,0.55) 
DepL[which(DepL==0)]=1 # Fixes "0" depths for Elkhorn Slough 
Kmax = rep(15,Npop); Kmax[5]=20
tmp = data.frame(Area = Areas,Estuary=SH4mn,Soft=SH1mn,Mixed=SH2mn,Hard=SH3mn,
                Depth=colMeans(DepH), Slope=DshHmn, Kelp=KlpHmn)
# Set params -----------------------------------------------------------------
# Set up Jags inputs -----------------------------------------------------
#
fitmodel = c("FitKmod.jags")
#  
jags.data <- list(Ncounts=Ncounts,NCarcsF=NCarcsF,NCarcsM=NCarcsM,
                  Psrv=Psrv,Ysrv=Ysrv,DemComp=DemComp,Npop=Npop,Nyrs=Nyrs,
                  AreaHD=AreaHD,AreaLD=AreaLD,NgrdH=NgrdH,NgrdL=NgrdL,
                  GridcountH=GridcountH,GridcountL=GridcountL,UB=UB,
                  DepH=DepH,DepL=DepL,KlpH=KlpH,KlpL=KlpL,propfem=propfem,
                  DshH=DshH,DshL=DshL,NoEstu=NoEstu,NoEstmn=NoEstmn,
                  SH1=SH1,SH2=SH2,SH3=SH3,SH4=SH4,HG=HG,LG=LG,
                  SL1=SL1,SL2=SL2,SL3=SL3,seKH=seKH,seKL=seKL,
                  SH1mn=SH1mn,SH2mn=SH2mn,SH3mn=SH3mn,SH4mn=SH4mn,
                  SL1mn=SL1mn,SL2mn=SL2mn,SL3mn=SL3mn,
                  KlpHmn=KlpHmn,DshHmn=DshHmn,KlpLmn=KlpLmn,DshLmn=DshLmn,
                  CarcShrkF=CarcShrkF,CarcShrkM=CarcShrkM,
                  PshkM=PshkM,YshkM=YshkM,PshkF=PshkF,YshkF=YshkF,
                  TotCarcF=TotCarcF,TotCarcM=TotCarcM,
                  LogInd=LogInd,LogPup=LogPup,
                  Ndeps=Ndeps,RelDensDep=RelDensDep,logN0=logN0,
                  PrpnDH=PrpnDH,PrpnDL=PrpnDL,Kmax=Kmax,
                  rho=rho,zeta=zeta,disprobF=disprobF,disprobM=disprobM) #CountP=CountP,CountI=CountI,  
#
inits <- function() list(sigT=runif(1, 1.5, 3), sigO=runif(1, .1, .2), 
                         sigK=runif(1, .3, .5), Nu = runif(1,.2,.5) )
#
params <- c("sigT","sigO","sigK","Nu","B0","B1","B2","B3","B4",
            "FRDF","UBE","Mbias","K","KHD","epsK","AdHz","N") 
#
nsamples <- 100
nthin <- 1
nadapt <- 100
nburnin <- 100
cores = detectCores()
ncore = min(20,cores-1)
#cl <- makeCluster(ncore)
nc <- ncore
#
# Run JAGS to fit model---------------------------------------------
#
out <- run.jags(data = jags.data, 
                monitor = params, 
                model = fitmodel,
                adapt=nadapt,
                inits = inits,
                n.chains = nc, 
                thin = nthin, 
                sample = nsamples, 
                burnin = nburnin,
                method="parallel") #inits = inits, 
#
# Diagnostic plots -------------------------------------------------
#
# stopCluster(cl)
post = rbind(out$mcmc[[1]], out$mcmc[[2]])
for (i in 3:nc){
  post = rbind(post, out$mcmc[[i]])
}
vn = varnames(out[["mcmc"]])
Nsims = dim(post)[1]
xx = which(!startsWith(vn,"N[") | !startsWith(vn,"AdHz[") )
post = post[,xx]
#
plot(out, vars = "sigT", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "sigO", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "sigK", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "FRDF", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "B0", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "B1", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "B2", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "B3", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "B4", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "Mbias", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "Nu", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "KHD", plot.type = c("trace", "histogram"),
     layout = c(1,2))
#
sumstats = summary(out)
#
plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",7]")),4],
    type="l", xlab = "year",ylab = "Otter abundance", 
    main = "Abundance trends, Monterey Peninsula")
iii = which(dfC$Pop==7)
points(dfC$Year[iii],dfC$Otts[iii])
plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",5]")),4],
     type="l", xlab = "year",ylab = "Otter abundance", 
     main = "Abundance trends, Elkhorn Slough")
iii = which(dfC$Pop==5)
points(dfC$Year[iii],dfC$Otts[iii])
plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",16]")),4],
     type="l", xlab = "year",ylab = "Otter abundance", 
     main = "Abundance trends, Pt Arguello")
iii = which(dfC$Pop==16)
points(dfC$Year[iii],dfC$Otts[iii])

AdHzMn = matrix(0,nrow=Npop,ncol=Nyrs)
PopReverse = seq(Npop,1,by=-1)
for (i in 1:Npop){
  for (j in 1:Nyrs) {
    # TemporalMn[i,j] = mean(outdf[,which(vn==paste0('Temporal[',i,',',j,']'))])
    AdHzMn[PopReverse[i],j] = sumstats[which(vn==paste0('AdHz[',i,',',j,']')),4]
  }
}
AdHzMn = t(AdHzMn)
layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(5,1), heights=c(1,1))
par(mar = c(3,5,2.5,2))
image(seq(1,Nyrs),seq(1,Npop),AdHzMn, col = heat.colors(100),
      xlab="Year",ylab="", yaxt="n", xaxt="n", cex.lab=.8)
title(main = "Spatiotemporal Variation in Shark Bite Hazard", font.main = 4)
# add categorical labels to y-axis
axis(1, at=seq(1,Nyrs), labels=as.character(c(1983:2017)),
     las=HORIZONTAL<-1, cex.axis=0.6)
axis(2, at=seq(1,Npop), labels=as.character(seq(Npop,1,by=-1)),
     las=HORIZONTAL<-1, cex.axis=0.6)
title(ylab = "Population", line = 2.5, cex.lab=.8)

ColorLevels <- seq(min(AdHzMn), max(AdHzMn), length=length(heat.colors(100)))
par(mar = c(3,2.5,2.5,2))
image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col = heat.colors(100),
      xlab="",ylab="",
      xaxt="n", cex.axis=0.8, cex.lab=.8)
title(main = "Color Scale", font.main = 1)
title(xlab = "log(Shark Bite Hazards)", line = .5, cex.lab=.8)
layout(1)

save.image(file="./Data/FitK_Results_temp.rdata")
