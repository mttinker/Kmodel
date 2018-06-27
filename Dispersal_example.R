library(extraDistr)
# Dispersal calculated using matrix multiplication separately for each sex,
# Using dispersal matrix M with #rows/collumns equal to number sub-pops
n = matrix(c(50,50,0),ncol=1)
M = matrix(c(.95,.04,.01,.025,.95,.025,.01,.04,.95),ncol=3,byrow = FALSE)
Mr = M; MA = M
nf =  matrix(nrow=3,ncol=100); nf[,1] = n
MalArea = matrix(data = 0, nrow=3,ncol=100); MalArea[3,1:50] = 1
for (i in 2:100){
  # For each pop, determine if Male area: if so, then dispersal to/from that area by females
  # is highly limited (reduced by factor of 100 here, could be a fitted param in model)
  # once transition from male area to female area then mean dispersal rates go to normal
  for (r in 1:3){
    if (MalArea[r,i] == 1){
      MA[r,] = .01*M[r,]
    }else{
      MA[r,] = M[r,]
    }
  }
  # Use dirichlet distribution to get stochastic probs for dispersal matrix
  for (c in 1:3){
    Mr[,c] = rdirichlet(1, 100*MA[,c])
  }
  # 
  nf[,i] = Mr%*%nf[,i-1]
}
matplot(t(nf), type = "l", xlab = "Year",ylab = "# females per sub-pop")

