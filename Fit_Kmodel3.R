# Script to fit a state-space population model with density-dependence
# to various data sets from California, including annual census counts,
# average abundance at regularly spaced 100m grid cells along the coast,
# And stranding data on cause of death ("Shark" vs "other").
# Load libraries -------------------------------------------------------------
require(runjags)
require(parallel)
require(gtools)
require(lattice)
require(coda)
require(ggplot2)
require(dplyr)
require(reshape2)
# Load data ------------------------------------------------------------------
# source("Process_data.R") # NOTE: need to run Process_data script first
rm(list = c())
load("../Data/Data_for_Kmod3.rdata")
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
sumstats = summary(out)

save.image(file="../Results/FitK_Results_mod3_test.rdata")

