#time varying case:
#X: private insurance
#L2: screening of rectal STI (yes/no)
#L3: log transformed CD4 (maximal value of 10)

datagen <- function(N, lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3,
                    beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4)
{
  ids <- seq(1,N,1)
  U <- runif(N, 0, 1) 
  U2 <- runif(N, 0, 1) 
  L2 = rbinom(N, 1, 0.5); 
  L1 = rbinom(N, 1, 0.5); 
  
  centres = c(1,2,3)
  X = rbinom(N, 2, plogis(-1+2*L2+L1))
  X = X+1
  A = rbinom(N, 1, plogis(alpha0+alpha1*X+0*L2+0*X*L2))
  denom = exp(beta0+beta1*A+beta2*X+beta3*L2-L1+2*A*L1)#-L1+4*L1*L2)
  Y = (-log(1-U)/denom)^(1/2)*(1/(lambda)^(1/2))
  denom2 = exp(-3.5-A+L2+X-L1)
  C = (-log(1-U2)/denom2)^(1/2)*(1/(0.1)^(1/2))
  #status = ifelse(Y<=12, 1, 0); 
  #status = ifelse(Y>C & C<=12, NA, status); 
  
  temp_data = data.frame(cbind(ids, X, L1, L2, A,C,Y))
  return(temp_data)
}


