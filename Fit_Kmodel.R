# Script to fit a state-space population model with density-dependence
# to various data sets from California, including annual census counts,
# average abundance at regularly spaced 100m grid cells along the coast,
# And stranding data on cause of death ("Shark" vs "other").
# Load libraries -------------------------------------------------------------
require(readxl)
require(stats)
require(MASS)
require(fitdistrplus)
require(gtools)
require(lme4)
require(lattice)
require(coda)
require(boot)
require(rjags)
require(runjags)
require(parallel)
require(doParallel)
require(ggplot2)
# Load data ------------------------------------------------------------------
# source("Process_data.R") # NOTE: need to run Process_data script first
load("Data_for_Kmod.rdata")
propfem = c(0.1,0.3,0.6)
# Set params -----------------------------------------------------------------
# Set up Jags inputs -----------------------------------------------------
#
fitmodel = c("FitKmod.jags")
#  
jags.data <- list(Ncounts=Ncounts,NCarcsF=NCarcsF,NCarcsM=NCarcsM,Npop=Npop,Nyrs=Nyrs,
                  Psrv=Psrv,Ysrv=Ysrv,CountI=CountI,DemComp=DemComp,
                  Area=Areas,AreaHD=AreaHD,AreaLD=AreaLD,NgrdH=NgrdH,NgrdL=NgrdL,
                  GridcountH=GridcountH,GridcountL=GridcountL,UB=UB,
                  DepH=DepH,DepL=DepL,KlpH=KlpH,KlpL=KlpL,propfem=propfem,
                  DshH=DshH,DshL=DshL,SubH=SubH,SubL=SubL,
                  CarcShrkF=CarcShrkF,CarcShrkM=CarcShrkM,
                  PshkM=PshkM,YshkM=YshkM,PshkF=PshkF,YshkF=YshkF,
                  rho=rho,zeta=zeta,disprobF=disprobF,disprobM=disprobM) #CountP=CountP,  
#
inits <- function() list(sigT=runif(1, 1.5, 3), sigE=runif(1, .2, .5), 
                         sigK=runif(1, .2, .5), Nu = runif(3,.2,.5) )
#
params <- c("sigT","sigE","sigK","Nu","B0","B1","B2","B3","B4",
            "ModalD","MAeff","K","AdHz") 
#
nsamples <- 100
nt <- 1
nb <- 500
cores = detectCores()
ncore = min(20,cores-1)
cl <- makeCluster(ncore)
nc <- ncore
#
# Run JAGS to fit model---------------------------------------------
#
out <- run.jags(data = jags.data, 
                monitor = params, 
                model = fitmodel, 
                inits = inits,
                n.chains = nc, 
                thin = nt, 
                sample = nsamples, 
                burnin = nb,
                method="rjparallel", cl=cl) #inits = inits, 
#
# Diagnostic plots -------------------------------------------------
#
stopCluster(cl)
post = rbind(out$mcmc[[1]], out$mcmc[[2]])
for (i in 3:nc){
  post = rbind(post, out$mcmc[[i]])
}
sumstats = summary(out)
vn = row.names(sumstats)
xx = which(!startsWith(vn,"N["))
post = post[,xx]
#
plot(out, vars = "sigT", plot.type = c("trace", "histogram"), layout = c(1,2))