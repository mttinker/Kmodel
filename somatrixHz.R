somatrixHz <- function(PrpK,sigma,AdHzF,AdHzM){
  #SOMATRIX function to generate 2x2 projection matrix for sea otters
  #  (classes = female/male) using proportional hazards formulation
  #   with density dependence, potential additive hazards (dens-independ, sex-specific) 
  #   and environmental stochasticity
  #
  # EXPLANATION:
  # Uses proportional hazards formulation to estimate survival
  #   Sx = exp(-gamma)
  #   gamma = exp(z0*sum(zeta1 + zeta2... + AdHz + eps))
  #     z0 = intercept (.01), ie. if no hazards 99% survival
  #     zeta1 = shared hazards for females and males when food abundant
  #     zeta2 = additional hazards for males
  #     zeta3 = Density-dependent hazards (scaled by proportion of K)
  #     AdHz = user supplied additional additive hazards 
  #        (which can be = 0 if no additional hazards)
  #     and eps = stochasticity, eps~rnorm(1,0,sigma)
  #     (if user-supplied sigma=0 then no stochasticity)
  #   R = reproductive output (weaned juveniles by females divided by 2, 50% male/female), 
  #      corrresponds to 1/2*br*wr*pM (pM=proportion mature females, ~0.65), and 
  #      uses power function to account for non-linear density-dependent effects
  #      (parameters were scaled to match telmetry estimates and basic logistic growth model,
  #       ~ 0.26 at max and ~0.18 at K, corresponds to max wr of 0.83, wr at K = 0.57)
  # NOTE1: user-supplied sigma represents std dev in repro output and log hazards,
  #    and is ~3-4x larger than the associated sd of log(lambda) 
  # NOTE2: hazard ratios are useful for fitting AdHz to carcass data (sex-sepcific)...
  #    ie ratio of exp(AdHz)/(exp(AdHz)+exp(sum(zeta))) gives expected proportion of carcasses from AdHz
  # NOTE3: expected pup count (assuming 75% pups born & still dependent at spring survey):
  #    pups = n[1,t-1]*fs*R*1.5
  #
  # OUTPUT:
  # lambda = algebraic value of lambda
  # ssd = associated stable stage distribution  
  # M = 4 by 4 projection matrix
  # vrts = vector of vital rates
  #
  # Zeta parameters for proportional hazards function (log hazards)
  z0 = .01 # Intercept (nuiscence parameter)
  z1 = 1.3	# Shared hazards for all animals when density low, food abundant
  z2 = 0.5  # Additional hazards for males
  z3 = 1.5  # Added density-dependent hazards (scaled by proportion of K)
  # Compute stochastic effects for this year 
  eps = rnorm(1,0,sigma) # env. stochasticity, variance in log hazard ratios
  epsR = exp(-1*(eps))/exp(.5*sigma^2) # Bias-adjusted log-normal error for R
  # Calculate vital rates
  R = min(.3,((2^(-3*PrpK))/10 + 0.167)*epsR)
  fs = exp(-z0*exp(z1+z3*PrpK+AdHzF+eps))
  ms = exp(-z0*exp(z1+z2+z3*PrpK+AdHzM+eps))
  M = matrix(c(fs*(1+R),   0,
               fs*R,      ms),byrow = TRUE,ncol = 2)
  lambda=eigen(M)$values[1]    # lambda=max value from vector of eigenvalues
  W=eigen(M)$vectors          # W=matrix of right eigenvectors 
  w=abs(W[,1])					      # w=stable distribution, unscaled
  ssd = w/sum(w)                # w=stable distribution, scaled
  result <- list(lam=lambda,SSD = ssd,M=M,
                 vrts = c(R,fs,ms))
  return(result)
}