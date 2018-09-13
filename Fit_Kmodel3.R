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
require(ggplot2)
require(dplyr)
require(reshape2)
# Load data ------------------------------------------------------------------
# source("Process_data.R") # NOTE: need to run Process_data script first
rm(list = c())
load("Data_for_Kmod3.rdata")
tauKg=1/seKg^2
Depthvect = seq(1,60); NDeps =  length(Depthvect)
# Set params -----------------------------------------------------------------
# Set up Jags inputs -----------------------------------------------------
#
fitmodel = c("FitKmod3.jags")
#  
jags.data <- list(Ncounts=Ncounts,NCarcsF=NCarcsF,Npop=Npop,Nyrs=Nyrs,
                  CountI=CountI,Sumcount=Sumcount,CarcShrkF=CarcShrkF,
                  TotCarcF=TotCarcF,PshkF=PshkF,YshkF=YshkF,Nbin=Nbin,
                  Pupratio=Pupratio,UB=UB,AreaB=AreaB,Areas=AreaHD,
                  Dep=Dep,Phrd=Phrd,Pklp=Pklp,Dsh=Dsh,Depthvect=Depthvect,
                  logN0=logN0,NDeps=NDeps,Psrv=Psrv,Ysrv=Ysrv)
#                  Kguess=Kguess,NKguess=NKguess,tauKg=tauKg 
#
inits <- function() list(sigT=runif(1, .45, .55), sigE=runif(1, .1, .3),
                         sigK=runif(1, .8, 1.2), Nu = runif(1,.5,1))
# sigK=runif(1, .1, .5),
params <- c("sigT","sigE","sigK","Nu","B1","B2","B3","B4","B5",
            "Intcpt","IntcptE","ModalD","Vinflate","Ntot",
            "ReproP","Sharkeff","UBE",
            "K","Kdns","AdHz","N") # 
#
nsamples <- 1000
nthin <- 10
nadapt <- 500
nburnin <- 2000
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
mcmc <- as.mcmc.list(out)
vn = varnames(out[["mcmc"]])
Nsims = length(as.matrix(mcmc)[,1])
Nsimsamp = sample(Nsims,1000,replace = TRUE)
#
plot(out, vars = "sig", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "Vinflate", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "Nu", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "ModalD", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "Intcpt", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "B", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "ReproP", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "Sharkeff", plot.type = c("trace", "histogram"),
     layout = c(1,2))
plot(out, vars = "UBE", plot.type = c("trace", "histogram"),
     layout = c(1,2))
# plot(out, vars = "K", plot.type = c("trace", "histogram"),
#      layout = c(1,2))
#
sumstats = summary(out)
#
plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",7]")),4],
    type="l", xlab = "year",ylab = "Otter abundance", 
    main = "Abundance trends, Monterey Peninsula",ylim = c(200,550))
iii = which(dfC$Pop==7)
points(dfC$Year[iii],dfC$Otts[iii])

plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",5]")),4],
     type="l", xlab = "year",ylab = "Otter abundance", 
     main = "Abundance trends, Elkhorn Slough",ylim = c(0,140))
iii = which(dfC$Pop==5)
points(dfC$Year[iii],dfC$Otts[iii])

plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",12]")),4],
     type="l", xlab = "year",ylab = "Otter abundance", 
     main = "Abundance trends, Cambria",ylim = c(200,500))
iii = which(dfC$Pop==12)
points(dfC$Year[iii],dfC$Otts[iii])

plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",13]")),4],
     type="l", xlab = "year",ylab = "Otter abundance", 
     main = "Abundance trends, Estero Bay",ylim=c(20,180))
iii = which(dfC$Pop==13)
points(dfC$Year[iii],dfC$Otts[iii])

plot(c(1983:2017),sumstats[which(startsWith(vn,"N[") & endsWith(vn,",15]")),4],
     type="l", xlab = "year",ylab = "Otter abundance", 
     main = "Abundance trends, Pismo Beach",ylim = c(0,250))
iii = which(dfC$Pop==15)
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

tempmat = as.matrix(mcmc)[Nsimsamp,startsWith(vn,"Kdns[")]
tempdf = as.data.frame(tempmat)
dfKdns = melt(tempdf)
dd= numeric()
for (i in 1:Npop){
#  ii = which(endsWith(as.character(dfKdns$variable),paste0("[",as.character(i),"]")))
  dd = which(endsWith(as.character(dfKdns$variable),paste0("[",as.character(i),"]")) & 
    dfKdns$value>1.1*sumstats[which(vn==paste0("Kdns[",as.character(i),"]")),3])
  dfKdns <- dfKdns[-dd,]
}
rm(tempmat,tempdf)
plt1 = ggplot(dfKdns, aes(x = variable, y = value)) +
  geom_violin(fill = "light blue", colour = "black",
              alpha = 0.7,scale = "width",trim = TRUE) +
  # geom_boxplot(fill = "light blue", colour = "black",
  #              alpha = 0.7) +
  scale_x_discrete(name = "Coastal Section", 
                   labels = as.character(seq(1,Npop))) +
  scale_y_continuous(name = "Density at K (otters per km2)",
                     limits=c(0,35)) +
  ggtitle("Posterior distributions, K density")
print(plt1)

Modal_Dep = as.matrix(mcmc)[Nsimsamp,which(vn=="ModalD")]
tempmat2 = as.matrix(mcmc)[Nsimsamp,which(startsWith(vn,"B"))]
tempmat = cbind(Modal_Dep,tempmat2)
tempdf = as.data.frame(tempmat)
dfbeta = melt(tempdf)
rm(tempmat,tempdf,Modal_Dep,tempmat2)
plt2 = ggplot(dfbeta, aes(x = variable, y = value)) +
  geom_violin(fill = "light blue", colour = "black",
              alpha = 0.7,scale = "width",trim = TRUE) +
  # geom_boxplot(fill = "light blue", colour = "black",
  #              alpha = 0.7) +
  scale_y_continuous(name = "Parameter Value") +
  scale_x_discrete(name = "Habitat Effect Parameter", labels = c("Modal Depth",
          "Depth Insh.","Depth Offsh.","Hard Substrate","Kelp Canopy","Shallow Slope")) +
#  ylim(-2,10) +
  ggtitle("Posterior distributions, Habitat Effect Parameters")
print(plt2)

CountI_Tot = numeric()
EstN_tot = numeric(); Est_LO = numeric(); Est_HI = numeric()
for (y in 1:Nyrs){
  CountI_Tot[y] = sum(dfC$Otts[dfC$yearN==y])
  if(CountI_Tot[y]==0){CountI_Tot[y]=NA}
  EstN_tot[y] = sumstats[vn = paste0("Ntot[",y,"]") ,4]
  Est_LO[y] = sumstats[vn = paste0("Ntot[",y,"]") ,1]
  Est_HI[y] = sumstats[vn = paste0("Ntot[",y,"]") ,3]
}
dfTrend = data.frame(Year = Years, Count = CountI_Tot,
                     Est_N = EstN_tot, LO = Est_LO, HI = Est_HI)
plt3 = ggplot(dfTrend) + 
  geom_line(aes(y=Est_N, x=Year, colour = "Estimate"),size=1)+
  geom_ribbon(aes(ymin=LO, ymax=HI, x=Year, fill = "95%CI"), alpha = 0.2)+
  geom_point(aes(x=Year,y=Count,colour="Survey Count"),size=2) +
  scale_colour_manual(name = '',  values=c("Estimate" = "blue", "Survey Count" = "black")) +
  scale_fill_manual(name = '',  values=c("95%CI" = "blue")) +
  labs(title = "Sea Otter Trends, State-Space Model Estimate vs Counts",
       x = "Year",y="Number Independent Otters")
print(plt3)

RelDens = numeric()
Est_LO = numeric(); Est_HI = numeric()
ModD = dfbeta$value[dfbeta$variable=="Modal_Dep"]
B1 = dfbeta$value[dfbeta$variable=="B1"]
B2 = dfbeta$value[dfbeta$variable=="B2"]
for (d in 1:max(Depthvect)){
  tempest = exp(-1*B1*0.01*pmax(0,ModD-Depthvect[d])^2 +
        -1*B2*0.01*max(0,Depthvect[d]-ModD)^2)
  RelDens[d] = mean(tempest)
  Est_LO[d] = quantile(tempest,c(.05))
  Est_HI[d] = quantile(tempest,c(.95))
}
scaleDep = 1/max(RelDens)
rm(tempest)
dfDepUse = data.frame(Depth = Depthvect,Relative_Density = RelDens*scaleDep,
                      LO = Est_LO*scaleDep, HI = Est_HI*scaleDep)
plt4 = ggplot(dfDepUse) + 
  geom_line(aes(y=Relative_Density, x=Depth, colour = "Mean"),size=1)+
  geom_ribbon(aes(ymin=LO, ymax=HI, x=Depth, fill = "95%CI"), alpha = 0.2)+
  scale_colour_manual(name = '',  values=c("Mean" = "blue")) +
  scale_fill_manual(name = '',  values=c("95%CI" = "blue")) +
  labs(title = "Sea Otter Relative Density by Depth",
       x = "Depth (m)",y="Relative Density at K")
print(plt4)

save.image(file="./Data/FitK_Results_mod3_test.rdata")

