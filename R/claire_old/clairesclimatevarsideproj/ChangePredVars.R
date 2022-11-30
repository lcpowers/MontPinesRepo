###### This code is trying to make it possible to run the Mont Pines jags model in a way that makes it easy to include/exclude particular predictor (particularly climate) variables.######

##### This should be able to happen by automatically generating inits where the inits for the predictor variables that should be excluded are essentially zero
##### Work flow: Prep data in 2LMER script, then run jags from this script

# Number of chains
nchains = 2


# Try bringing prior for clouds here
# This does not work
# cloud_dist = paste("grwth_cloudCoef ~ dnorm(0,1e-10)")

# This does not work. Tried with underscore and . in the name. Does not work when you add the variable name to the data call at the end of the jags model. Tried putting cloud.dist in eval() & in eval(as.name())
cldmin = -1e10
cldmax = 1e10

fogmin = -1e10
fogmax = 1e10

tmayseptmin = -1e10
tmayseptmax = 1e10

toctaprmin = -1e10
toctaprmax = 1e10

precipmin = -1e10
precipmax = 1e10

# # Inits maybe don't matter too much, but this seems to work
# TempMaySept_init <- c(0.01,0.001) # 1e-10 ## Ask about what distribution would make sense for these
# TempOctApr_init <- c(0.01,0.001)
# Cloud_init <- rep(1e-10,nchains)
# Fog_init <- c(0.01,0.001)
# PrecipAugJuly_init <- c(0.01,0.001)

jags.mod <- run.jags("R/NoRepro_JAGSmodel_Feb1.R", # Call to specific jags model
                     n.chains=2,
                     data=mp4jags,
                     burnin=500,
                     thin=10,
                     sample=3000,
                     adapt=500,
                     method="parallel")

jags.mod.ChangePreVars.NoCloud <- jags.mod


### Notes:
# - Model runs when define init value for a particular climate variable before inits, then use that variable to assign values to inits
#     - Did not check whether output was reasonable
# - Model also seems to run when defined climate vars in this script right before running
# - Model runs when climate vars defined above, and when a single value is selected from a distribution for inits that we want to be non-zero

vec <- rexp(n=10000,rate = 3)
hist(vec)


