## Doak Lab Monterey Pines Project
## Data preparation and call to JAGS script
## Modify data and assign variables needed for the JAGS model with data lags (missing years)
## Info about associated JAGS script later
##
#### This file is the equivalent of `erbr_2lmerStepsToJAGS_210504.R`
#### I'm going through that file line-by-line, bringing over the code that seems relevant to us. 
##

## Note: Details on mixed models given here: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
## Other relevant links:
# https://bayesball.github.io/BOOK/bayesian-multiple-regression-and-logistic-models.html#bayesian-logistic-regression
# http://www.columbia.edu/~cjd11/charles_dimaggio/DIRE/styled-4/styled-11/code-8/#logistic-regression

rm(list=ls())
# graphics.off()

## LOAD PACKAGES AND FUNCTION --------------------------------------------------------------------
library(lme4)
library(tidyverse)
library(ggplot2)
library(rjags)
library(runjags)
library(dplyr)
library(coda)
library(corrplot)
## SET WD ON LOCAL COMPUTER (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------
setwd("./R/") # should work if within project but will give error if already in 'R' folder or elsewhere

## LOAD DATA --------------------------------------------------------------------------------------
montpines <- read.csv("../data/fullannual data C.csv") %>% 
  select(-1)  ## Get rid of first column that seems like old row numbers from elsewhere

# montpines <- montpines[1:763,]
# head(montpines)

montpines$growth <- NA # Size this year minus size last measured year - skips lag>1 years
montpines$propgrowth <- NA # (size this year minus size last year)/size last measured year skips lag >1 years
montpines$allgrowth <- NA # Growth for all years -- deals with lag >1 years
montpines$allpropgrowth <- NA # Proportional growth for all years -- deals with lag >1 years

montpines$sitetag <- paste0(montpines$Site,montpines$TagNo)
montpines <- arrange(montpines,sitetag,demoyr)
montpines$yr <- montpines$demoyr

for(i in 2:nrow(montpines)) {
  
  # Total growth and proportional growth, not dividing across years without data
  if(montpines$sitetag[i]==montpines$sitetag[i-montpines$lags[i]]){ # if this row and lag row have the same sitetag (are the same individual)
    
    ### Growth not accounting for/filling in years without data ###
    # Total growth = current DBH - last measured DBH
    montpines$growth[i] = montpines$DBH[i] - montpines$DBH[i-montpines$lags[i]] 
    
    # Proportional growth = total growth/last measured DBH
    montpines$propgrowth[i] = montpines$growth[i]/montpines$DBH[i-montpines$lags[i]] 
   
  if(!is.na(montpines$surv[i]) & montpines$surv[i]==1 & montpines$lags[i]==0) {
    montpines$growth[i] = NA
    montpines$propgrowth[i] = NA
  }
  
    ### Growth accounting for/filling in years without data ###
  # If there is data for two consecutive years (lag=1), don't divide growth
  if(montpines$lags[i]==1) {
    montpines$allgrowth[i] = montpines$growth[i]
    montpines$allpropgrowth[i] = montpines$propgrowth[i]
  }
  # If the lag is greater than 1, distribute growth between all skipped years
  if(montpines$lags[i]>1){ 

    l = montpines$lags[i]
    growthvec = rep(montpines$growth[i]/l,l)
    montpines$allgrowth[(i-l+1):i] = growthvec
    propgrowthvec = rep(montpines$propgrowth[i]/l,l)
    montpines$allpropgrowth[(i-l+1):i] = propgrowthvec

    }
  }
  
}

rm(propgrowthvec,growthvec,i,l)

## SETTING UP NECESSARY VARIABLES -----------------------------------------------------------------
##### Setting up the jags model with lagged values

# number of rows in data
Nallrows <- length(montpines$Site)
numyears <- n_distinct(montpines$yr)
numsites <- n_distinct(montpines$Site)

# Create goodgrowrows from new lag str size variable. 
# lagsrtsz should be equal to lags except in the case where 
# plt died and the lag is > 0. In this case lagsrtsz = 0. 
montpines$lagsrtsz = NA
for (i in 1:nrow(montpines)) { 
  if (montpines[i,]$lags > 0 & montpines[i,]$surv == 0){ 
    montpines$lagsrtsz[i] = 0
  }
  else {
    lagval = montpines[i,]$lags
    montpines$lagsrtsz[i] = lagval}
}

## Identify rows that are good dependent values (ending sizes) for surv or growth  
goodrows <- which(montpines$lags>0) # This finds the rows with data
goodgrowrows <- which(montpines$lagsrtsz > 0) # This finds rows with good growth data

lagvals <- montpines$lags[goodrows]
Ncases <- length(goodrows)
Ngrowcases <- length(goodgrowrows)
Survs <- montpines$surv
dbh <- montpines$DBH


##### Variable setup for repro fitting #####
rows.w.sz <- which(!is.na(montpines$DBH)&montpines$resurveyarea==1) # rows with size values that are also in resurv area
rows.wo.sz <- which(is.na(montpines$DBH)&montpines$resurveyarea==1) # rows without size values also in resurv area
Ndirectszcases <- length(rows.w.sz)
Nindirectszcases <- length(rows.wo.sz)
rows.w.recruits <- which(montpines$newplt > 0)

## create vector to indicate if alive or dead after missing yr(s)
montpines$RowNum <- 1:nrow(montpines) # Add a column to indicate row number
rows.wo.sz.alive <- as.data.frame(matrix(NA, # initiate data.frame
                                         nrow=length(rows.wo.sz),
                                         ncol=2))
colnames(rows.wo.sz.alive) <- c("Rows", "Alive") 
rows.wo.sz.alive$Rows <- rows.wo.sz # insert row indexes for rows without size

for (i in rows.wo.sz) {                                                  # Loop over all tags with 1 or more yrs of missing dbh data               
  tag.val <- montpines$sitetag[i]                                          # get tag no 
  tag.each <- subset(montpines, montpines$sitetag==tag.val)                # Subset main data for tag i
  tag.surv <- tag.each$surv[!is.na(tag.each$surv) & tag.each$RowNum>i]   # Get non-missing years of survival after first year
  rows.wo.sz.alive$Alive[rows.wo.sz.alive$Rows==i] <- tag.surv[1]        # Store surv for 1st non-missing yr post each missed yr
  rm(tag.val,tag.each,tag.surv)
}

rows.wo.sz.alive$Alive[is.na(rows.wo.sz.alive$Alive)] <- 0  # Change NAs to 0, these are lines where missed yr was last & recorded as dead
rows.wo.sz.alive <- rows.wo.sz.alive$Alive                  # Change to vector

## At this point in eriogonum code there are lines to make transects random effects. Skipping these for now based on our conversation where Dan said it probably made more sense to treat sites as fixed effects (and no transects in the mont pines data)

## Make yr and transect numerical to use in jags as random effects 
montpines$Year.num <- as.factor(montpines$demoyr) %>% as.numeric()
Year.num <- montpines$Year.num

montpines$transect.num <- as.factor(montpines$Site) %>% as.numeric()
transect.num <- montpines$transect.num

# linear index of transect-year combos 
yrtranscombo=100*montpines$transect.num+montpines$Year.num 

## Set up a logical variable that is whether there is surv or grwth data in a yr (0,1): this is needed for the summing up of infl nums to predict new plts
# datayesno <- rep(1,length(montpines$surv))
# datayesno[which(is.na(montpines$surv)==TRUE)] <- 0 

## ------------------------------------------------------------------------------------------------

## Make dataframe w new plts (that are likely recent seedlings) for each transect & yr ------------
## Make df that will hold data containing new plants

years <- unique(montpines$Year.num)
years <- years[order(years)]
montpines.newPlts <- as.data.frame(rep(unique(montpines$transect.num), each=length(years)))
colnames(montpines.newPlts) <- "TransectNew.num" 
montpines.newPlts$Year.num <- rep(years)

newPlts <- montpines[montpines$newplt == 1, ] # Identify new plants 
num.newPlts <- newPlts %>%                    # Count N of newplts per yr&trans
     group_by(transect.num,yr) %>% 
     summarise(num.newPlts=n())
colnames(num.newPlts) = c("TransectNew.num", "Year.num", "num.newPlts") #to match April's code 

## Add number of new plants to df of each transect & year
montpines.newPlts <- left_join(montpines.newPlts, num.newPlts, by=c("TransectNew.num", "Year.num"))
montpines.newPlts$num.newPlts[is.na(montpines.newPlts$num.newPlts)] <- 0   #Change NAs (no new plants) to zeros
montpines.newPlts$num.newPlts[montpines.newPlts$Year.num==1] <- NA         #Change new plts in year 1 to NA

## Add column so new plts in t+1 match year t. Note: variables w '1' at the end represent t+1 data
montpines.newPlts <- montpines.newPlts %>% mutate(num.newPlts1=lead(num.newPlts))  

## Add climate variables to new plants data 
montpines.climate <- montpines %>% select(c(Year.num, Temp.Oct.Apr., 
                                            Temp.May.Sept., Precip.Aug.July, fog, cloud)) 
clim <- unique(montpines.climate)
montpines.newPlts <- left_join(montpines.newPlts, clim, by="Year.num")
#should this line (below) be lead or lag? (April had this as 'lead,' ours should be lag, right?)
montpines.newPlts <- montpines.newPlts %>% mutate(TempOctApr1 = lead(Temp.Oct.Apr.), 
                                                  TempMaySept1 = lead(Temp.May.Sept.), 
                                                  PrecipAugJuly1 = lead(Precip.Aug.July),
                                                  fog1 = lead(fog), 
                                                  cloud1 = lead(cloud)) 

## Remove 2 lines that correspond to transect-year combos that are not in the main data file
montpines.newPlts$yrtranscombo=100*montpines.newPlts$TransectNew.num+montpines.newPlts$Year.num
yrtrans.unq <- unique(yrtranscombo)
montpines.newPlts <- montpines.newPlts[montpines.newPlts$yrtranscombo %in% yrtrans.unq,]

montpines.newPlts <- montpines.newPlts[is.na(montpines.newPlts$num.newPlts1) == FALSE,] #I am not sure why this is necessary but in erio code...
newplts <- montpines.newPlts$num.newPlts1
newplt.trans <- montpines.newPlts$TransectNew.num
newplt.yr <- montpines.newPlts$Year.num
newPltlines <- length(montpines.newPlts$TransectNew.num)

## Make a linear index of transect-year combos for new plts
newplt.yrtranscombo=100*newplt.trans+newplt.yr 


## do we need to add repro variables to this mp4jags? Nope
mp4jags <- montpines %>% 
  select(demoyr,DBH,surv,
         lags,lagsrtsz,
         TempOctApr=Temp.Oct.Apr.,
         TempMaySept=Temp.May.Sept.,
         PrecipAugJuly=Precip.Aug.July,
         fog,cloud,
         transect.num)

###########################################

## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------

modname <- ""
jags.mod <- run.jags(paste0('montpines_JAGSmodel',modname,'.R'), # Call to specific jags model
                     n.chains=3,
                     data=mp4jags,
                     burnin=500,
                     thin=10,
                     sample=3000,
                     adapt=5000,
                     method="parallel")

failed.jags('model')


save(jags.mod, file=paste0('modeloutput/mp_JAGSmod_',modname,"_",Sys.Date(),'.rdata'))
saveRDS(jags.mod, paste0('modeloutput/mp_JAGSmod_',modname,"_",Sys.Date(),'.rds'))


## ------------------------------------------------------------------------------------------------

# 
# ## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
# jags.mod <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_210613.rds")
summary(jags.mod)
plot(jags.mod)
summ.mod <- summary(jags.mod)
tail(summ.mod[,1:5], n=32)
gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)

chains <- jags.mod$mcmc
chains <- bind_rows(lapply(chains, as.data.frame))
colMeds <- apply(chains,2,median)
colSDs <- apply(chains,2,sd)


## Look at correlation b/w params
chains.1 <- chains %>% dplyr::select(!contains(c("randomeffect", "precision")))
chains.1 <- chains.1 %>% dplyr::select(!c(deviance, resid.sum.sq))

cor.chains <- cor(chains.1)
corrplot(cor.chains, method="circle", type="lower")


## ** Make bar graph comparing median param ests & 80% CIs b/w diff datasets **
## ------------------------------------------------------------------------------------------------


