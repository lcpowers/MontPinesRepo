######################################################################################################
######################################################################################################
#### JAGS model file written by runjags version 2.0.4-6 on 2020-10-23 09:28:28 
######################################################################################################
######################################################################################################

## RosNew == DBH

### Model template as follows - ensure this is syntactically correct before running the model!

model{

# In the BUGS/JAGS language we must use an explicit for loop:

#########################################################################
## Doing one loop to generate the estimated sz for each plt in each yr w a survival data point (all these have sz, except the surv=0 yrs) 
## This loop also estimates the survival to each yr: for unobserved yrs, this includes the chances of surv each unobserved yr in the past 
## Cases in this loop are pairs of observation yrs, that can include yrs without any observations 
for(i in 1:Ncases){

  ## Do first lag, which must use the starting size (DBH)
  regression_mean[goodrows[i]-lagvals[i]+1] <- grwth_intercept 
                                              + grwth_dbhCoef*DBH[goodrows[i]-lagvals[i]]                  # growth coefficient
                                               + grwth_TempMaySeptCoef*TempMaySept[goodrows[i]-lagvals[i]]     # Summer T coef
     #                                          + grwth_TempOctAprCoef*TempOctApr[goodrows[i]-lagvals[i]]       # Winter T coef 
                                               + grwth_PrecipAugJulyCoef*PrecipAugJuly[goodrows[i]-lagvals[i]] # Annual Precip coef
     #                                          + grwth_cloudCoef*cloud[goodrows[i]-lagvals[i]]                 # cloud coef
                                               + grwth_fogCoef*fog[goodrows[i]-lagvals[i]]                     # fog coef
      + grwth_Transect_randomeffect[transect.num[goodrows[i]-lagvals[i]]]
      
  r.growth[goodrows[i]-lagvals[i]+1] <- exp(grwthvar_intercept + grwthvar_dbhCoef*DBH[goodrows[i]-lagvals[i]]) 
  
  ## Survival prob, based on the last years observed size 
  Surv_mu[goodrows[i]-lagvals[i]+1] <- 1/(1+exp(-(surv_intercept
                                                  + surv_dbhCoef*DBH[goodrows[i]-lagvals[i]]
                                                   + surv_TempMaySeptCoef*TempMaySept[goodrows[i]-lagvals[i]]
                                                  # + surv_TempOctAprCoef*TempOctApr[goodrows[i]-lagvals[i]]
                                                   + surv_PrecipAugJulyCoef*PrecipAugJuly[goodrows[i]-lagvals[i]]
                                                  # + surv_cloudCoef*cloud[goodrows[i]-lagvals[i]]
                                                   + surv_fogCoef*fog[goodrows[i]-lagvals[i]]
                                                   + surv_Transect_randomeffect[transect.num[goodrows[i]-lagvals[i]]]))) 
                                                  
     
  ## This is the loop of remaining lagged yrs for this row i: note that the : operator is diff in jags than r & can't be decreasing, hence the use of negative lag: also, if lag=1, then this will be from 0 to -1, & the loop will be skipped 
  ## Survival rates are based on the inferred previous year size, since it was not observed 
  for (j in (goodrows[i]-lagvals[i]+1):(goodrows[i]-1)) { # For the consecutive rows of no data that correspond to good row i # tried removing +1 and -1 from range
    
     regression_mean[j+1] <- grwth_intercept 
                             + grwth_dbhCoef*regression_mean[j]      # growth coefficient
                              + grwth_TempMaySeptCoef*TempMaySept[j]     #  Summer T coef
                             # + grwth_TempOctAprCoef*TempOctApr[j]       # Winter T coef 
                              + grwth_PrecipAugJulyCoef*PrecipAugJuly[j] # Annual Precip coef
                             # + grwth_cloudCoef*cloud[j]                 # cloud coef
                              + grwth_fogCoef*fog[j]                     # fog coef
                              + grwth_Transect_randomeffect[transect.num[j]] # transect(site) random effect
                           
     r.growth[j+1] <- exp(grwthvar_intercept + grwthvar_dbhCoef*regression_mean[j])

     Surv_mu[j+1] <- Surv_mu[j]*1/(1+exp(-(surv_intercept 
                                           + surv_dbhCoef*regression_mean[j]
                                            + surv_TempMaySeptCoef*TempMaySept[j]
                                           # + surv_TempOctAprCoef*TempOctApr[j]
                                           # + surv_PrecipAugJulyCoef*PrecipAugJuly[j]
                                           # + surv_cloudCoef*cloud[j]
                                            + surv_fogCoef*fog[j]
                                            + surv_Transect_randomeffect[transect.num[j]])))
                                            
   } #End lags loop 
} #End of going through cases for growth and survival
#####################################################################

  #####################################################################
  ## A loop that makes each repro amount estimate for years with observed sizes. 
  ## maxRepro is the max repro output of largest individuals 
  for(i in 1:Ndirectszcases){
    
   
    #repro_size[rows.w.sz[i]] <- repro_intercept #intercept only model 
    
    repro_size[rows.w.sz[i]] <- max(0, reproDBHcoef*DBH[rows.w.sz[i]] + repro_intercept)
    # repro_size[rows.w.sz[i]] <- maxRepro/(1+exp(-(reproDBHcoef*DBH[rows.w.sz[i]] + repro_intercept))) #  +
                                                    # repro_TempMaySeptCoef*TempMaySept1[rows.w.sz[i]] + 
                                                    # repro_TempOctAprCoef*TempOctApr1[rows.w.sz[i]] + 
                                                    # repro_PrecipAugJulyCoef*PrecipAugJuly1[rows.w.sz[i]] + 
                                                    # repro_cloudCoef*cloud1[rows.w.sz[i]] + 
                                                    # repro_fogCoef*fog1[rows.w.sz[i]] + 
                                                    # repro_Transect_randomeffect[transect.num[rows.w.sz[i]]])))  
  }
  
  # repro_size[rows.w.sz[i]] <- max(0, reproDBHcoef*DBH[rows.w.sz[i]] + repro_intercept)
  ############t#########################################################
  
  #####################################################################
  # Do the same for repro estimates, for years without an observation, so based on inferred size
  for(i in 1:Nindirectszcases){
    
    #repro_size[rows.w.sz[i]] <- repro_intercept*Surv_mu #intercept only model 
    
   repro_size[rows.wo.sz.alive[i]] <- Surv_mu[rows.wo.sz.alive[i]]*max(0, reproDBHcoef*regression_mean[rows.wo.sz.alive[i]] + repro_intercept)
   #  repro_size[rows.wo.sz.alive[i]] <- Surv_mu[rows.wo.sz.alive[i]]*maxRepro/(1+exp(-(reproDBHcoef*regression_mean[rows.wo.sz.alive[i]] + repro_intercept))) # +
                                                           # repro_TempMaySeptCoef*TempMaySept1[rows.wo.sz.alive[i]] +
                                                           # repro_TempOctAprCoef*TempOctApr1[rows.wo.sz.alive[i]] +
                                                           # repro_PrecipAugJulyCoef*PrecipAugJuly1[rows.wo.sz.alive[i]] +
                                                           # repro_cloudCoef*cloud1[rows.wo.sz.alive[i]] +
                                                           # repro_fogCoef*fog1[rows.wo.sz.alive[i]] +
                                                           # repro_Transect_randomeffect[transect.num[rows.wo.sz.alive[i]]])))
  }
  
  # repro_size[rows.wo.sz.alive[i]] <- Surv_mu[rows.wo.sz.alive[i]]*max(0, reproDBHcoef*regression_mean[rows.wo.sz.alive[i]] + repro_intercept)
  #####################################################################
  #this is a little loop to maybe help debug repro fitting - it is throwing an error with some undefined repro rows ? 
  for (i in 1:Nexcludedcases){
    repro_size[rows.excluded.repro[i]] <- 0
  }
  
  #####################################################################
  ## This loop is getting the predicted num of tot recruits for each transect/year combo, 
  # & fitting to num of new plts a yr later. 
  # This should be sum of reprosizes for that trans*yr 
  # The inprod & the 1st logical are making 0,1 for whether a given line is part of each trans-yr combo
  
  for (i in 1:newPltlines){
    
    pred.tot.recruits[i] <- inprod(yrtranscombo==newplt.yrtranscombo[i], repro_size) 
    #r.newplts[i] <- pred.tot.recruits[i] #is this right? we just want the mean prediction? 
    
    #p.newplts[i] <- r.newplts/(r.newplts+pred.tot.recruits[i]) #think this was for neg bin? 
    # newplts[i] ~  dnegbin(p.newplts[i], r.newplts)
    newplts[i] ~ dpois(pred.tot.recruits[i]) 
  }
  #####################################################################
  
  
  
  #####################################################################
  ## Then, a loop that matches the estimated sizes with the observed sizes, just for rows with observed sizes
  for(i in 1:Ngrowcases){
    ## These lines describe the response distribution and linear model terms
    
    regression_residual[i] <- DBH[goodgrowrows[i]] - regression_mean[goodgrowrows[i]]
    DBH[goodgrowrows[i]] ~ dnorm(regression_mean[goodgrowrows[i]], 1/r.growth[goodgrowrows[i]])
    
  } #End of going through cases for sizes
  #####################################################################
  
  #####################################################################
  ## Now a loop to do survival predictions, for the ending years with observations
  for(i in 1:Ncases){
    ## Fit the survival function:
    Survs[goodrows[i]] ~ dbern(Surv_mu[goodrows[i]])
  } #End loop to do survival predictions
  #####################################################################

## do we need to fit estimated repro to actual repro anywhere???
  
########################################    
## These lines give the prior distributions for the parameters to be estimated
# r_disp ~ dgamma(0.01, 0.01) # this is the dispersion for all cases, as is standard
r_disp ~ dunif(0,50)
# Note: the prior for the dispersion parameter k is quite important for convergence
# [A DuMouchel prior may be better than a Gamma prior]

## GROWTH  
grwth_intercept ~ dunif(-5,5) ## Do these values make sense for DBH?
grwth_dbhCoef ~ dunif(-3,3) ## Do these values make sense for DBH?
grwth_TempMaySeptCoef ~ dunif(tmayseptmin,tmayseptmax) ## Do these values make sense for DBH?
grwth_TempOctAprCoef ~ dunif(toctaprmin,toctaprmax) ## Do these values make sense for DBH?
grwth_PrecipAugJulyCoef ~ dunif(precipmin,precipmax) ## Do these values make sense for DBH?
grwth_fogCoef ~ dunif(fogmin,fogmax) ## Do these values make sense for DBH?
grwth_cloudCoef ~ dunif(cldmin,cldmax) ## Do these values make sense for DBH? 

# ## Growth transect random effects
for(transect.num_iterator in 1:5){
  grwth_Transect_randomeffect[transect.num_iterator] ~ dnorm(0, grwth_Transect_precision) }
grwth_Transect_precision ~ dunif(1,10)


## VARIANCE IN GROWTH
grwthvar_intercept ~ dunif(-3,3) 
grwthvar_dbhCoef ~ dunif(-3,3) 


## SURVIVAL
surv_intercept ~  dnorm(0, 10^-6) 
surv_dbhCoef ~  dnorm(0, 10^-6) 

surv_TempMaySeptCoef ~ dnorm(0,tmayseptmax)
# surv_TempMaySeptCoef ~ dnorm(0, 10^-6) 

surv_TempOctAprCoef ~ dnorm(0, toctaprmax)
# surv_TempOctAprCoef ~ dnorm(0, 10^-6) 

surv_PrecipAugJulyCoef ~ dnorm(0,precipmax)
# surv_PrecipAugJulyCoef ~ dnorm(0, 10^-6) 

surv_fogCoef ~ dnorm(0, fogmax) 
# surv_fogCoef ~ dnorm(0, 10^-6) 

surv_cloudCoef ~ dnorm(0, cldmax) 
# surv_cloudCoef ~ dnorm(0, 10^-6) 

## Survival transect random effects
for(transect.num_iterator in 1:5){
  surv_Transect_randomeffect[transect.num_iterator] ~ dnorm(0, surv_Transect_precision) }
surv_Transect_precision ~ dunif(0,1)

## REPRODUCTION
maxRepro ~ dunif(0, 100) 
reproDBHcoef ~ dnorm(0, 10^-6) 
repro_intercept ~ dnorm(0, 10^-6) 
r.newplts ~ dnorm(0.01,0.01)
# repro_size ~ dunif(0,10)

repro_TempMaySeptCoef ~ dnorm(0, tmayseptmax)
# repro_TempMaySeptCoef ~ dnorm(0, 10^-6)

repro_TempOctAprCoef ~ dnorm(0, toctaprmax)
# repro_TempOctAprCoef ~ dnorm(0, 10^-6)

repro_PrecipAugJulyCoef ~ dnorm(0, precipmax)
# repro_PrecipAugJulyCoef ~ dnorm(0, 10^-6)

repro_fogCoef ~ dnorm(0, fogmax)
# repro_fogCoef ~ dnorm(0, 10^-6)

repro_cloudCoef ~ dnorm(0, cldmax)
# repro_cloudCoef ~ dnorm(0, 10^-6)

## Repro transect random effects 
for(transect.num_iterator in 1:5){
  repro_Transect_randomeffect[transect.num_iterator] ~ dnorm(0, repro_Transect_precision) }
repro_Transect_precision ~ dunif(0,1)

#resid.sum.sq <- sum(regression_residual^2)
} # end of model specification 

## do we need to add anything to monitor?

## ADD BACK TO MONITOR (optional): dic
# These lines are hooks to be read by runjags (they are ignored by JAGS):
#monitor# deviance,grwth_intercept,grwth_Transect_randomeffect,surv_Transect_randomeffect,repro_Transect_randomeffect,grwth_dbhCoef,grwth_TempMaySeptCoef,grwth_TempOctAprCoef,grwth_PrecipAugJulyCoef,grwth_fogCoef,grwth_cloudCoef,grwthvar_intercept,grwthvar_dbhCoef,surv_intercept,surv_dbhCoef,surv_TempMaySeptCoef,surv_TempOctAprCoef,surv_PrecipAugJulyCoef,surv_fogCoef,surv_cloudCoef,maxRepro,reproDBHcoef,repro_intercept,repro_TempMaySeptCoef,repro_TempOctAprCoef,repro_PrecipAugJulyCoef,repro_fogCoef,repro_cloudCoef
#modules# glm on
#response# DBH
#residual# regression_residual
#fitted# regression_fitted
#data# Survs,Ncases,Ngrowcases,Ndirectszcases,Nindirectszcases,Nexcludedcases,goodrows,goodgrowrows,lagvals,DBH,TempMaySept,TempOctApr,PrecipAugJuly,fog,cloud,transect.num,newPltlines,rows.wo.sz.alive,rows.w.sz,rows.excluded.repro,yrtranscombo,newplt.yrtranscombo,TempOctApr1,TempMaySept1,PrecipAugJuly1,fog1,cloud1,cldmin,cldmax,fogmin,fogmax,precipmin,precipmax,tmayseptmin,tmayseptmax,toctaprmin,toctaprmax

######################################################################################################
#### Initial values 
######################################################################################################
######################################################################################################

inits{
  "grwth_intercept" <- 1/10
  "grwth_dbhCoef" <- 1/10
  "grwth_TempMaySeptCoef" <- 0
  "grwth_TempOctAprCoef" <- 0
  "grwth_PrecipAugJulyCoef" <- 0
  "grwth_fogCoef" <- 0
  "grwth_cloudCoef" <-0
  "grwth_Transect_precision" <- 2
  
  "grwthvar_intercept" <- 0.01
  "grwthvar_dbhCoef" <- 0.01
  
  "surv_intercept" <- -1/10
  "surv_dbhCoef" <- 1/10
  "surv_TempMaySeptCoef" <- 0
  "surv_TempOctAprCoef" <- 0
  "surv_PrecipAugJulyCoef" <- 0
  "surv_fogCoef" <- 0
  "surv_cloudCoef" <- 0
  "surv_Transect_precision" <- 0.001
  
  "maxRepro" <- 0.01 
  "repro_intercept" <- 0.01
  "reproDBHcoef" < -1/10
  "repro_TempMaySeptCoef" <- 0
  "repro_TempOctAprCoef" <- 0
  "repro_PrecipAugJulyCoef" <- 0
  "repro_fogCoef" <- 0
  "repro_cloudCoef" <- 0
  "repro_Transect_precision" <- 0.01
  "r.newplts" <- 1
  # "repro_size" <- 0.01
}

# inits{
# 
#   "grwth_intercept" <- 1.5
#   "grwth_dbhCoef" <- 1
#   "grwth_TempMaySeptCoef" <- 0
#   "grwth_TempOctAprCoef" <- 0
#   "grwth_PrecipAugJulyCoef" <- 0
#   "grwth_fogCoef" <- 0
#   "grwth_cloudCoef" <- 0
#   "grwth_Transect_precision" <- 2
# 
#   "grwthvar_intercept" <- 0.1
#   "grwthvar_dbhCoef" <- 0.1
# 
#   "surv_intercept" <- -10
#   "surv_dbhCoef" <- -0.01
#   "surv_TempMaySeptCoef" <- 0
#   "surv_TempOctAprCoef" <- 0
#   "surv_PrecipAugJulyCoef" <- 0
#   "surv_fogCoef" <- 0
#   "surv_cloudCoef" <- 0
#   "surv_Transect_precision" <- 0.001
#   
#   "maxRepro" <- 0.01 
#   "repro_intercept" <- 0.01
#   "reproDBHcoef" < -1
#   "repro_TempMaySeptCoef" <- 0
#   "repro_TempOctAprCoef" <- 0
#   "repro_PrecipAugJulyCoef" <- 0
#   "repro_fogCoef" <- 0
#   "repro_cloudCoef" <- 0
#   "repro_Transect_precision" <- 0.01
#   "r.newplts" <- 1
#   # "repro_size" <- 0.01
# 
# }
# inits{
# 
#   "grwth_intercept" <- 1.5
#   "grwth_dbhCoef" <- 1
#   "grwth_TempMaySeptCoef" <- 0
#   "grwth_TempOctAprCoef" <- 0
#   "grwth_PrecipAugJulyCoef" <- 0
#   "grwth_fogCoef" <- 0
#   "grwth_cloudCoef" <- 0
#   "grwth_Transect_precision" <- 2
# 
#   "grwthvar_intercept" <- 0.01
#   "grwthvar_dbhCoef" <- -0.01
# 
#   "surv_intercept" <- -20
#   "surv_dbhCoef" <- 0.01
#   "surv_TempMaySeptCoef" <- 0
#   "surv_TempOctAprCoef" <- 0
#   "surv_PrecipAugJulyCoef" <- 0
#   "surv_fogCoef" <- 0
#   "surv_cloudCoef" <- 0
#   "surv_Transect_precision" <- 0.001
# 
#   "maxRepro" <- 0.01 #this seems unreasonable 
#   "repro_intercept" <- 0.01
#   "reproDBHcoef" < -1
#   "repro_TempMaySeptCoef" <- 0
#   "repro_TempOctAprCoef" <- 0
#   "repro_PrecipAugJulyCoef" <- 0
#   "repro_fogCoef" <- 0
#   "repro_cloudCoef" <- 0
#   "repro_Transect_precision" <- 0.01
#   "r.newplts" <- 1
#   # "repro_size" <- 0.01
# }