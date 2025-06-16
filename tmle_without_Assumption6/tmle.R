
## CREATE DATASET MYDATA.CSV
  library(geepack);library(MASS);library(ResourceSelection);library(ltmle); library(SuperLearner); library(survival); library(coxphw); library(rootSolve);
  library(dplyr); library(glm2); library(nleqslv);
  library(data.table); library(splines); library(nnet);
  setDTthreads(1)
  #library(reshape2)  #do not use for data frame only
  
  logit <- function(term) {
    return( ifelse(!is.na(term),log(term/(1-term)),NA) )
  }
  
  EXPIT <- function(term) {
    return( ifelse(!is.na(term),exp(term)/(1+exp(term)),NA) )
  }
  
  source("datagen.R")
  set.seed(100001)
  seeds = floor(runif(1000)*10^8);
  set.seed(seeds[1]) ## SEED USED TO GENERATE THE mydata.csv DATAFRAME
  
  N <- 2500
  alpha0=-1; alpha1=1; alpha2=0; alpha3=0; 
  lambda = 5
  beta0=-4; beta1=-2; beta2=1; beta3=-1;
  df <- datagen(N, lambda, alpha0=alpha0, alpha1=alpha1, alpha2=alpha2, alpha3=alpha3,
                beta0=beta0, beta1=beta1, beta2=beta2, beta3=beta3, beta4=beta4)
  
  
  #### READ IN DATA
  df = read.csv("mydata.csv", header=T)
  tau = 5
  df$C = ifelse(df$C<=tau+1, df$C, tau+1) #administrative censoring at time 6
  txt = 1 ### HERE IS WHERE YOU INDICATE TREATMENT OR CONTROL
  sourcenumber = 1 ### WHICH DATA SOURCE WE ARE LOOKING AT
  
  ### Treatment model and Centre probability
  afit = glm2(A ~ X, family = binomial(link=logit), data = df)  ## HERE YOU SPECIFY TREATMENT MODEL
  df$pred_obs = predict(afit, newdata = df, type="response")
  df$pred_obs = ifelse(df$A==1, df$pred_obs, 1-df$pred_obs)

  ## first estimate the hazards with cox regression
  df$T_tilde = pmin(df$Y, df$C)
  df$delta = ifelse(df$T_tilde==df$Y & df$Y<=tau, 1, 0)
  df$AL1 = df$A*df$L1
  fit_init = coxph(Surv(T_tilde, delta) ~ A + L1 + L2 + X + AL1, data = df, method="breslow") ## HERE YOU SPECIFY YOUR HAZARDS MODEL FOR EVENT OF INTEREST
  baseline <- basehaz(fit_init,centered=F)
  baseline$start = baseline$time
  baseline$end = c(baseline$time[2:nrow(baseline)], tau+1+0.00001)
  coef_init = fit_init$coefficients
  df$z = as.vector(as.matrix(df %>% select(A,L1,L2,X,AL1)) %*% coef_init) ## ESTIMATED PARAMETRIC COMPONENT IN h(t|Z), i.e., exp(Z'beta)
  df$surv_tau = exp(-baseline[baseline$start<=tau & baseline$end>tau,]$hazard*exp(df$z))
  
  ## censoring
  df$C_stat = ifelse(df$T_tilde==df$C & df$C<=tau, 1, 0)
  fit_c = coxph(Surv(T_tilde, C_stat) ~ A + X + L2 + L1, data = df, method="breslow") ## HERE YOU SPECIFY YOUR CENSORING HAZARDS MODEL
  coef_c = fit_c$coefficients
  df$z_c = as.matrix(df %>% select(A,X,L2,L1)) %*% coef_c   ## ESTIMATED PARAMETRIC COMPONENT IN h_c(t|Z)
  baselinec <- basehaz(fit_c,centered=F)
  baselinec$start = baselinec$time
  baselinec$end = c(baselinec$time[2:nrow(baselinec)], tau+1+0.00001)
  
  tmpdat = df
  tmpdat = tmpdat[order(tmpdat$T_tilde),]
  score = NULL
  
  ## CREATE ANOTHER DATASET WHERE TREATMENT IS SET TO TXT TO PREDICT FINAL SURVIVAL PROBABILITY WITH A=txt
  impdat = tmpdat
  impdat$A=txt; impdat$AL1 = impdat$A*impdat$L1;
  impdat$z = as.vector(as.matrix(impdat %>% select(A,L1,L2,X,AL1)) %*% coef_init) # SAME COVARATES AS IN "FIT_INIT" ABOVE, BUT WITH A FIXED AT A=txt
  tmpdat$surv_tau_a = exp(-baseline[baseline$start<=tau & baseline$end>tau,]$hazard*exp(impdat$z))
  
  #########################################
  #########################################
  ### DO NOT CHANGE CODE STARTING FROM HERE
  #########################################
  #########################################

  dH = NULL; start=NULL; # here I calculate the baseline hazards
  for(i in 1:nrow(tmpdat))
  {
    tmptime = tmpdat$T_tilde[i]
    Expected_num_i = ifelse(tmpdat$T_tilde==tmptime & tmpdat$delta==1, 1, 0)
    Expected_denom_i = ifelse(tmpdat$T_tilde>=tmptime, exp(tmpdat$z),0)
    Expected_num = sum(Expected_num_i)
    Expected_denom = sum(Expected_denom_i)
    dHtmp = Expected_num/Expected_denom
    dH  = c(dH, dHtmp)
    start = c(start, tmptime)
  }
  
  hazard = cumsum(dH) #H_0(t)
  baseline = as.data.frame(cbind(hazard, start))
  baseline$end = c(baseline$start[2:nrow(baseline)], tau+1+0.0001)
  baseline$dH = dH
  
  surv_c <- lapply(tmpdat$T_tilde, function(t) {
    exp(-baselinec[baselinec$start<=t & baselinec$end>t,]$hazard*exp(tmpdat$z_c))} )
  
  myfunc = function(gamma)
  {
    score=0
    for(i in 1:nrow(tmpdat))
    {
      tmptime = tmpdat$T_tilde[i]
      tmpdat$surv_c = surv_c[[i]]
      tmpdat$surv_t = exp(-baseline[baseline$start<=tmptime & baseline$end>tmptime,]$hazard*exp(tmpdat$z))
      tmpdat$weight_out = I(tmpdat$T_tilde>=tmptime & tmpdat$A==txt & tmpdat$X==sourcenumber)*(tmpdat$surv_tau)/(tmpdat$pred_obs*tmpdat$surv_c*tmpdat$surv_t)
      
      tmpdat$weight_in = (tmpdat$surv_tau)/(tmpdat$surv_t)
      tmpdat$hazt = baseline[baseline$start<=tmptime & baseline$end>tmptime,]$dH*exp(tmpdat$z)*exp(tmpdat$weight_in*gamma)
      
      tmpdat$stat = 0; tmpdat$stat[i] = ifelse(tmpdat$delta[i]==1, 1, 0);
      blah = tmpdat$weight_out*(tmpdat$stat - tmpdat$hazt)
      score_i = ifelse(tmptime<=tau, sum(blah), 0); 
      score  = c(score, score_i)
    }
    sum(score)
  } 
  mygamma = 0
  gamma = nleqslv(mygamma, myfunc, control=list(allowSingular=FALSE))$x
  
  ## reiterate #calculate gamma by writing own code based on algorithm in paper
  tmp = which(baseline$start<=tau & baseline$end>tau)
  myfuncup = function(newgamma)
  {
    score=0
    for(i in 1:nrow(tmpdat))
    {
      tmptime = tmpdat$T_tilde[i]
      tmpdat$surv_c = surv_c[[i]]
      tmpdat$surv_t = exp(-cumhaz_1[[i]])
      
      tmpdat$weight_out = I(tmpdat$T_tilde>=tmptime & tmpdat$A==txt & tmpdat$X==sourcenumber)*(tmpdat$surv_tau_up)/(tmpdat$pred_obs*tmpdat$surv_c*tmpdat$surv_t)
      tmpdat$weight_in = tmpdat$surv_tau_up/tmpdat$surv_t
      tmpdat$hazt = dH1[[i]]*exp(tmpdat$weight_in*newgamma)
      
      tmpdat$stat = 0; tmpdat$stat[i] = ifelse(tmpdat$delta[i]==1, 1, 0);
      blah = tmpdat$weight_out*(tmpdat$stat - tmpdat$hazt)
      score_i = ifelse(tmptime<=tau, sum(blah), 0); 
      score  = c(score, score_i)
    }
    sum(score)
  }
  
  cumhaz_1_a <- list() ## cumhaz(1) builds on weight_0
  for (llen in 1:nrow(tmpdat)) {
    cumhaz_1_a[[llen]] <- baseline[baseline$start<=tmpdat$T_tilde[llen] & baseline$end>tmpdat$T_tilde[llen],]$hazard*exp(impdat$z)
  }
  
  count=0
  while(abs(gamma)>1e-03 & count<=30)
  {
    ### update
    if(count==0){
      weight_0 = list() ## create weights for time 0 used for dH1
      for(i in 1:nrow(tmpdat))
      {
        tmptime = tmpdat$T_tilde[i]
        tmpdat$surv_t = exp(-baseline[baseline$start<=tmptime & baseline$end>tmptime,]$hazard*exp(tmpdat$z))
        weight_0[[i]] = tmpdat$surv_tau/tmpdat$surv_t
      }
      dH1 <- list() ## dH(1) builds on weight_0
      for (llen in 1:nrow(tmpdat)) {
        dH1[[llen]] <- baseline[baseline$start<=tmpdat$T_tilde[llen] & baseline$end>tmpdat$T_tilde[llen],]$dH*exp(tmpdat$z)*exp(weight_0[[llen]]*gamma) 
      }
      cumhaz_1 <- Reduce(`+`, dH1, accumulate = TRUE) ## cumhaz(1) builds on weight_0
      weight_0_a = list() ## create weights for time 0 used for dH1
      for(i in 1:nrow(tmpdat))
      {
        tmptime = tmpdat$T_tilde[i]
        tmpdat$surv_t = exp(-baseline[baseline$start<=tmptime & baseline$end>tmptime,]$hazard*exp(impdat$z))
        weight_0_a[[i]] = tmpdat$surv_tau_a/tmpdat$surv_t
      }
      dH1_a <- list() ## dH(1) builds on weight_0
      for (llen in 1:nrow(tmpdat)) {
        dH1_a[[llen]] <- baseline[baseline$start<=tmpdat$T_tilde[llen] & baseline$end>tmpdat$T_tilde[llen],]$dH*exp(impdat$z)*exp(weight_0_a[[llen]]*gamma) 
      }
      cumhaz_1_a <- Reduce(`+`, dH1_a, accumulate = TRUE)
      tmpdat$surv_tau_up = exp(-cumhaz_1[[tmp]]) ## surv_tau (1)
      tmpdat$surv_tau_up_a = exp(-cumhaz_1_a[[tmp]])
    } else{
      weight_up=list()
      for(i in 1:nrow(tmpdat))
      {
        tmptime = tmpdat$T_tilde[i]
        tmpdat$surv_t = exp(-cumhaz_1[[i]])
        weight_up[[i]] = tmpdat$surv_tau_up/tmpdat$surv_t
      }
      for (llen in 1:nrow(tmpdat)) {
        dH1[[llen]] <- dH1[[llen]]*exp(weight_up[[llen]]*gamma) }
      cumhaz_1 <- Reduce(`+`, dH1, accumulate = TRUE)
      weight_up_a=list()
      for(i in 1:nrow(tmpdat))
      {
        tmptime = tmpdat$T_tilde[i]
        tmpdat$surv_t = exp(-cumhaz_1_a[[i]])
        weight_up_a[[i]] = tmpdat$surv_tau_up_a/tmpdat$surv_t
      }
      for (llen in 1:nrow(tmpdat)) {
        dH1_a[[llen]] <- dH1_a[[llen]]*exp(weight_up_a[[llen]]*gamma) }
      cumhaz_1_a <- Reduce(`+`, dH1_a, accumulate = TRUE)
      tmpdat$surv_tau_up = exp(-cumhaz_1[[tmp]])
      tmpdat$surv_tau_up_a = exp(-cumhaz_1_a[[tmp]])
    }
    
    ## reiterate
    #calculate gamma by writing own code based on algorithm in paper
    mygamma=0;
    gammaup = nleqslv(mygamma, myfuncup, control=list(allowSingular=F))$x
    gamma=gammaup
    gamma
    count=count+1
  }
  
  ## outcome model results
  tmpdat$surv_tau_final = exp(-cumhaz_1_a[[tmp]])
  fulldat = tmpdat[tmpdat$X==sourcenumber,]
  mean(fulldat$surv_tau_final)
  
  #final estimate
  tk_double = mean(fulldat$surv_tau_final)

