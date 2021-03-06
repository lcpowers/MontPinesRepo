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
                                              + grwth_TempOctAprCoef*TempOctApr[goodrows[i]-lagvals[i]]       # Winter T coef 
                                              + grwth_PrecipAugJulyCoef*PrecipAugJuly[goodrows[i]-lagvals[i]] # Annual Precip coef
                                              + grwth_cloudCoef*cloud[goodrows[i]-lagvals[i]]                 # cloud coef
                                              + grwth_fogCoef*fog[goodrows[i]-lagvals[i]]                     # fog coef
     + grwth_Transect_randomeffect[transect.num[goodrows[i]-lagvals[i]]]
      
  r.growth[goodrows[i]-lagvals[i]+1] <- exp(grwthvar_intercept + grwthvar_dbhCoef*DBH[goodrows[i]-lagvals[i]]) 
  
  ## Survival prob, based on the last years observed size 
  Surv_mu[goodrows[i]-lagvals[i]+1] <- 1/(1+exp(-(surv_intercept
                                                  + surv_dbhCoef*DBH[goodrows[i]-lagvals[i]]
                                                  + surv_TempMaySeptCoef*TempMaySept[goodrows[i]-lagvals[i]]
                                                  + surv_TempOctAprCoef*TempOctApr[goodrows[i]-lagvals[i]]
                                                  + surv_PrecipAugJulyCoef*PrecipAugJuly[goodrows[i]-lagvals[i]]
                                                  + surv_cloudCoef*cloud[goodrows[i]-lagvals[i]]
                                                  + surv_fogCoef*fog[goodrows[i]-lagvals[i]]
                                                  + surv_Transect_randomeffect[transect.num[goodrows[i]-lagvals[i]]]))) ### What do we want to do with ? ###
                                                  
     
  ## This is the loop of remaining lagged yrs for this row i: note that the : operator is diff in jags than r & can't be decreasing, hence the use of negative lag: also, if lag=1, then this will be from 0 to -1, & the loop will be skipped 
  ## Survival rates are based on the inferred previous year size, since it was not observed 
  for (j in (goodrows[i]-lagvals[i]+1):(goodrows[i]-1)) { # For the consecutive rows of no data that correspond to good row i # tried removing +1 and -1 from range
    
     regression_mean[j+1] <- grwth_intercept 
                             + grwth_dbhCoef*regression_mean[j]         # growth coefficient
                             + grwth_TempMaySeptCoef*TempMaySept[j]     #  Summer T coef
                             + grwth_TempOctAprCoef*TempOctApr[j]       # Winter T coef 
                             + grwth_PrecipAugJulyCoef*PrecipAugJuly[j] # Annual Precip coef
                             + grwth_cloudCoef*cloud[j]                 # cloud coef
                             + grwth_fogCoef*fog[j]                     # fog coef
                             + grwth_Transect_randomeffect[transect.num[j]] # transect(site) random effect
                           
     r.growth[j+1] <- exp(grwthvar_intercept + grwthvar_dbhCoef*regression_mean[j])

     Surv_mu[j+1] <- Surv_mu[j]*1/(1+exp(-(surv_intercept 
                                           + surv_dbhCoef*regression_mean[j]
                                           + surv_TempMaySeptCoef*TempMaySept[j]
                                           + surv_TempOctAprCoef*TempOctApr[j]
                                           + surv_PrecipAugJulyCoef*PrecipAugJuly[j]
                                           + surv_cloudCoef*cloud[j]
                                           + surv_fogCoef*fog[j]
                                           + surv_Transect_randomeffect[transect.num[j]])))
                                          
   } #End lags loop 
} #End of going through cases for growth and survival
#####################################################################

  #####################################################################
  ## A loop that makes each reproyesno and repro amount estimate for years with observed sizes. 
  ## maxRepro is the max repro output of largest individuals 
  for(i in 1:Ndirectszcases){
  
    #Dan has a = maxRepro, b = bRepro, c = cRepro - I just renamed so we can keep track more easily..?
    repro_size[rows.w.sz[i]] <- maxRepro/(1+exp(-(bRepro+cRepro*DBH[i])))
  }
  #####################################################################
  
  #####################################################################
  ## Do the same for repro estimates, for years without an observation, so based on inferred size
  for(i in 1:Nindirectszcases){
    
    # repro_size is same as loop with direct size cases, except we need to multiply by 
    # probability of being alive. is regression_mean what we want here instead of DBH? 
    # I am not sure that indexing Surv_mu and regression_mean will give us what we want..? 
    # do we want rows.wo.sz.alive or rows.wo.size? (since we have surv_mu in there???) 
    repro_size[rows.wo.sz.alive[i]] <- maxRepro/(1+exp(-(bRepro+cRepro*regression_mean[i]))) * 
                                          Surv_mu[i]
    
  }
  #####################################################################
  
  ## Loop over only rows with >0 recruits, to compare estimated to observed. 
  # Do we need this loop?
  for(i in 1:Nrows.w.recruits) {
    
    p.recruits[i] <- r.recruits/(r.recruits + repro_size[rows.w.recruits[i]])
    newPltlines[rows.w.recruits[i]] ~ dnegbin(p.recruits[i], r.recruits)
 
  }
  
  
  #####################################################################
  ## This loop is getting the predicted num of tot recruits for each transect/year combo, 
  # & fitting to num of new plts a yr later. 
  # This should be sum of reprosizes for that trans*yr, * a term for climate effects in the year prior 
  # The inprod & the 1st logical are making 0,1 for whether a given line is part of each trans-yr combo
  
  for (i in 1:newPltlines){
    
    pred.tot.recruits[i] = inprod(yrtranscombo==newplt.yrtranscombo[i], reprosize) #is this getting a sum? 
    sum.repro.sizes = NULL 
    ### Ask dan about these lines . Do we need both new_plt intercept and pred.tot.recruits??? 
    #check lag vars 
    mn.new.plts[i] = exp(newplt_intercept + pred.tot.recruits[i]+
                           repro_TempMaySeptCoef*Temp_MaySept1[i] +
                           repro_TempOctAprCoef*TempOctApr1[i] + 
                           repro_PrecipAugJulyCoef*PrecipAugJuly1[i] + 
                           repro_cloudCoef*cloud1[i]+
                           repro_fogCoef*fog1[i]) #didn't include random effect in here yet 
    
    p.newplts[i] <- r.newplts/(r.newplts+mn.new.plts[i])
    newplts[i] ~  dnegbin(p.newplts[i], r.newplts)
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
  

########################################    
## These lines give the prior distributions for the parameters to be estimated
# r_disp ~ dgamma(0.01, 0.01) # this is the dispersion for all cases, as is standard
r_disp ~ dunif(0,50)
# Note: the prior for the dispersion parameter k is quite important for convergence
# [A DuMouchel prior may be better than a Gamma prior]

## GROWTH  
grwth_intercept ~ dunif(-5,5) 
grwth_dbhCoef ~ dunif(-3,3) 
grwth_TempMaySeptCoef ~ dunif(-2,2)
grwth_TempOctAprCoef ~ dunif(-2,2)
grwth_PrecipAugJulyCoef ~ dunif(-2,2)
grwth_fogCoef ~ dunif(-2,2)
####
grwth_cloudCoef ~ dunif(cldmin,cldmax)
####

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
surv_TempMaySeptCoef ~ dnorm(0, 10^-6)
surv_TempOctAprCoef ~ dnorm(0, 10^-6)
surv_PrecipAugJulyCoef ~ dnorm(0, 10^-6) 
surv_fogCoef ~ dnorm(0, 10^-6)
surv_cloudCoef ~ dnorm(0, 10^-6)

## Survival transect random effects
for(transect.num_iterator in 1:5){
  surv_Transect_randomeffect[transect.num_iterator] ~ dnorm(0, surv_Transect_precision) }
surv_Transect_precision ~ dunif(0,1)

## NEW PLANTS
maxRepro ~ dnorm(0, 10^-6)
bRepro ~ dnorm(0, 10^-6)
cRepro ~ dnorm(0, 10^-6)
r.recruits ~ dnorm(0, 10^-6)
# r.newplts ~ dgamma(0.01,0.01)


#resid.sum.sq <- sum(regression_residual^2)
} # end of model specification 



## ADD BACK TO MONITOR (optional): dic
# These lines are hooks to be read by runjags (they are ignored by JAGS):
#monitor# deviance,grwth_intercept,grwth_Transect_randomeffect,surv_Transect_randomeffect,grwth_dbhCoef,grwth_TempMaySeptCoef,grwth_TempOctAprCoef,grwth_PrecipAugJulyCoef,grwth_fogCoef,grwth_cloudCoef,grwthvar_intercept,surv_intercept,surv_dbhCoef,surv_TempMaySeptCoef,surv_TempOctAprCoef,surv_PrecipAugJulyCoef,surv_fogCoef,surv_cloudCoef
#modules# glm on
#response# DBH
#residual# regression_residual
#fitted# regression_fitted
#data# Ncases,Ngrowcases,goodrows,goodgrowrows,lagvals,DBH,TempMaySept,TempOctApr,PrecipAugJuly,fog,cloud,transect.num,cldmin,cldmax




######################################################################################################
#### Initial values 
######################################################################################################
######################################################################################################

inits{
  "grwth_intercept" <- 2
  "grwth_dbhCoef" <- 1
  "grwth_TempMaySeptCoef" <- 0.01
  "grwth_TempOctAprCoef" <- 0.01
  "grwth_PrecipAugJulyCoef" <- 0.01
  "grwth_fogCoef" <- 0.01
  "grwth_cloudCoef" <- 0.01
  "grwth_Transect_precision" <- 2
  
  "grwthvar_intercept" <- 0.01
  "grwthvar_dbhCoef" <- 0.01
  
  "surv_intercept" <- -1
  "surv_dbhCoef" <- 0.01
  "surv_TempMaySeptCoef" <- 0.01
  "surv_TempOctAprCoef" <- 0.01
  "surv_PrecipAugJulyCoef" <- 0.01
  "surv_fogCoef" <- 0.01
  "surv_cloudCoef" <- 0.01
  "surv_Transect_precision" <- 0.001
# 
  # "reproyesno_intercept" <- -1
  # "reproyesno_RosCoef" <- 0.1
  # "reproyesno_TempFallCoef" <- 0.01
  # "reproyesno_TempSummerCoef" <- 0.01
  # "reproyesno_TempWinterCoef" <- 0.01   
  # "reproyesno_PptFallCoef" <- 0.01
  # "reproyesno_PptSummerCoef" <- 0.01  
  # "reproyesno_Transect_precision" <- 0.01
 # "repro_precision" <- 0.01
#  "repro_intercept" <- 2
  "repro_RosCoef" <- 1
  "repro_TempMaySeptCoef" <- 0.01
  "repro_TempOctAprCoef" <- 0.01
  "repro_PrecipAugJulyCoef" <- 0.01
  "repro_fogCoef" <- 0.01
  "repro_cloudCoef" <- 0.01
 # "repro_Transect_precision" <- 0.01
 # 
 # "newplt_intercept" <- -1
 # "r.newplts" <- 1
}

inits{

  "grwth_intercept" <- 1.5
  "grwth_dbhCoef" <- 1
  "grwth_TempMaySeptCoef" <- 0.001
  "grwth_TempOctAprCoef" <- 0.001
  "grwth_PrecipAugJulyCoef" <- 0.001
  "grwth_fogCoef" <- 0.001
  "grwth_cloudCoef" <- 0.001
  "grwth_Transect_precision" <- 2

  "grwthvar_intercept" <- 0.1
  "grwthvar_dbhCoef" <- 0.1

  "surv_intercept" <- -10
  "surv_dbhCoef" <- -0.01
  "surv_TempMaySeptCoef" <- 0.05
  "surv_TempOctAprCoef" <- 0.5
  "surv_PrecipAugJulyCoef" <- 0.5
  "surv_fogCoef" <- 0.5
  "surv_cloudCoef" <- 0.5
  "surv_Transect_precision" <- 0.001

  #   "reproyesno_intercept" <- -1
  #"repro_RosCoef" <- 1
  "repro_TempMaySeptCoef" <- 0.01
  "repro_TempOctAprCoef" <- 0.01
  "repro_PrecipAugJulyCoef" <- 0.01
  "repro_fogCoef" <- 0.01
  "repro_cloudCoef" <- 0.01
  #   "reproyesno_Transect_precision" <- 0.01
  #
  #   "repro_precision" <- 0.01
  #   "repro_intercept" <- 2
  #   "repro_RosCoef" <- 1
  #   "repro_PptFallCoef" <- 0.01
  #   "repro_PptSummerCoef" <- 0.01
  #   "repro_TempWinterCoef" <- 0.01
  #   "repro_TempFallCoef" <- 0.01
  #   "repro_TempSummerCoef" <- 0.01
  #   "repro_Transect_precision" <- 0.01
  #
  #   "newplt_intercept" <- -0.1
  #   "r.newplts" <- 1
}
inits{

  "grwth_intercept" <- 1.5
  "grwth_dbhCoef" <- 1
  "grwth_TempMaySeptCoef" <- 0.001
  "grwth_TempOctAprCoef" <- 0.001
  "grwth_PrecipAugJulyCoef" <- 0.001
  "grwth_fogCoef" <- 0.001
  "grwth_cloudCoef" <- 0.001
  "grwth_Transect_precision" <- 2

  "grwthvar_intercept" <- 0.01
  "grwthvar_dbhCoef" <- -0.01

  "surv_intercept" <- -20
  "surv_dbhCoef" <- 0.01
  "surv_TempMaySeptCoef" <- 0.05
  "surv_TempOctAprCoef" <- 1
  "surv_PrecipAugJulyCoef" <- 0.5
  "surv_fogCoef" <- 0.5
  "surv_cloudCoef" <- 0.5
  "surv_Transect_precision" <- 0.001

  # "reproyesno_intercept" <- -5
  # "reproyesno_RosCoef" <- 0.1
  # "reproyesno_TempFallCoef" <- 0.01
  # "reproyesno_TempSummerCoef" <- 0.01
  # "reproyesno_TempWinterCoef" <- 0.01
  # "reproyesno_PptFallCoef" <- 0.01
  # "reproyesno_PptSummerCoef" <- 0.01
  # "reproyesno_Transect_precision" <- 0.01
  #
  # "repro_precision" <- 0.01
  # "repro_intercept" <- 2
  # "repro_RosCoef" <- 0.5
  # "repro_PptFallCoef" <- 0.01
  # "repro_PptSummerCoef" <- 0.01
  # "repro_TempWinterCoef" <- 0.01
  # "repro_TempFallCoef" <- 0.01
  # "repro_TempSummerCoef" <- 0.01
  # "repro_Transect_precision" <- 0.01
  #
  # "newplt_intercept" <- 0
  # "r.newplts" <- 1
}