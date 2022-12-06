## Doak Lab Monterey Pines Project
## Data preparation and call to JAGS script
## Modify data and assign variables needed for the JAGS model with data lags (missing years)##
#### This file is the equivalent of `erbr_2lmerStepsToJAGS_210504.R`

## Note: Details on mixed models given here: https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf
## Other relevant links:
# https://bayesball.github.io/BOOK/bayesian-multiple-regression-and-logistic-models.html#bayesian-logistic-regression
# http://www.columbia.edu/~cjd11/charles_dimaggio/DIRE/styled-4/styled-11/code-8/#logistic-regression

rm(list=ls())

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
if(basename(getwd())!="R") setwd("./R/") # should work if within project but will give error if already in 'R' folder or elsewhere

montpines <- read_csv("../data/fullannual data C.csv") %>% 
  select(-1) %>% ## Get rid of first column that seems like old row numbers from elsewhere
  # filter(demoyr < 2016) %>% 
  mutate(growth = NA, # Size this year - size last measured year - skips lag > 1 years
         propgrowth = NA, # (size this year - size last year)/size last measured year. Skips lag >1 years
         allgrowth = NA,  # Growth for all years -- deals with lag >1 years
         allpropgrowth = NA, # Proportional growth for all years -- deals with lag >1 years
         sitetag = paste0(Site,TagNo),
         yr = demoyr) %>%
  arrange(sitetag,yr)  %>% 
  select(-c(demoyr,DFN,dbh1,ht1,dbhfromDFN))

# Get climate vals per site per yr
climate <- montpines %>% 
  select(Site,yr,Temp.Oct.Apr.,Temp.May.Sept.,Precip.Aug.July,fog,cloud) %>% 
  unique() 

## Growth loop --------------------------------
# This loop finds growth from size measurments. Included is growth from one census to the next (timeframe varies) as well as annual growth (total growth between censuses/years btwn censuses) as well as proportional total and annual growth
for(i in 2:nrow(montpines)) {
  
  # Total growth and proportional growth, not dividing across years without data
  if(montpines$sitetag[i]==montpines$sitetag[i-montpines$lags[i]]){ # if this row and lag row have the same sitetag (are the same individual)...
    
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
    # If there is data for two consecutive years (lag==1), don't divide growth
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

# total number of rows in data
Nallrows <- nrow(montpines)
numyears <- n_distinct(montpines$yr)
numsites <- n_distinct(montpines$Site)

# Create goodgrowrows from new lag str size variable. 
# lagsrtsz should be equal to lags except in the case where 
# plt died and the lag is > 0. In this case lagsrtsz = 0. 
montpines$lagsrtsz = NA
for (i in 1:nrow(montpines)) { 
  
  if (montpines$lags[i] > 0 & montpines$surv[i] == 0){ # if good data in row i & tree did not survive then 0
    
    # Set the lag value equal to zero
    montpines$lagsrtsz[i] = 0
    
  }
  else {
    
    # Else, keep lag value the same as in lags col
    montpines$lagsrtsz[i] = montpines$lags[i]}
}

## Identify rows that are good dependent values (ending sizes) for surv or growth  
goodrows <- which(montpines$lags>0) # This finds the rows with data (size and survival outcome)
goodgrowrows <- which(montpines$lagsrtsz > 0) # This finds rows with good growth data (size when survival outcome == 1)

lagvals <- montpines$lags[goodrows] # Lag values where there is data
Ncases <- length(goodrows) # Number of rows with good data
Ngrowcases <- length(goodgrowrows) # Number of rows with good growth data
Survs <- montpines$surv # Isolate survival column?
dbh <- montpines$DBH # Isolate DBH column?

# Compute lagged values for environmental/climate vars. 
# This climate cols so that last years climate vars are associated with this yr's size 
montpines$TempMaySept1 = NA
montpines$TempOctApr1 = NA
montpines$PrecipAugJuly1 = NA
montpines$fog1 = NA
montpines$cloud1 = NA

for(i in 1:nrow(montpines)){
  
  yr = montpines$yr[i]
  site = montpines$Site[i]
  
  # This skips 2002 for now, until we (maybe) have 2001 clim data
  montpines$TempMaySept1[i] = ifelse((yr-1)%in%unique(montpines$yr[montpines$Site==site]), # If climate data for last year exists for this site
                                     climate$Temp.May.Sept.[climate$Site==site&climate$yr==(yr-1)], # Use that value for the lagged clim val
                                     NA) # else, use NA
  
  # Same logic as above for all other clim vars
  montpines$TempOctApr1[i] = ifelse((yr-1)%in%unique(montpines$yr[montpines$Site==site]),
                                    climate$Temp.Oct.Apr.[climate$Site==site&climate$yr==(yr-1)],NA)
  
  montpines$PrecipAugJuly1[i] = ifelse((yr-1)%in%unique(montpines$yr[montpines$Site==site]),
                                       climate$Precip.Aug.July[climate$Site==site&climate$yr==(yr-1)],NA)
  
  montpines$fog1[i] = ifelse((yr-1)%in%unique(montpines$yr[montpines$Site==site]),
                             climate$fog[climate$Site==site&climate$yr==(yr-1)],NA)
  
  montpines$cloud1[i] = ifelse((yr-1)%in%unique(montpines$yr[montpines$Site==site]),
                               climate$cloud[climate$Site==site&climate$yr==(yr-1)],NA)
  
}

##### Variable setup for repro fitting #####
rows.w.sz <- which(!is.na(montpines$DBH)&montpines$resurveyarea==1) # rows with size values that are also in resurv area
rows.wo.sz <- which(is.na(montpines$DBH)&montpines$resurveyarea==1) # rows without size values also in resurv area
Ndirectszcases <- length(rows.w.sz) # number of rows with measured sizes & in area surveyed for recruits
rows.w.recruits <- which(montpines$newplt==1) # recruit rows
#####

## create vector to indicate if alive or dead after missing yr(s)
montpines$RowNum <- 1:nrow(montpines) # Add a column to indicate row number
rows.wo.sz.alive <- as.data.frame(matrix(NA, # initiate data.frame
                                         nrow=length(rows.wo.sz),
                                         ncol=2))
colnames(rows.wo.sz.alive) <- c("Rows", "Alive")
rows.wo.sz.alive$Rows <- rows.wo.sz # insert row indexes for rows without size

for (i in rows.wo.sz) {                                                  # Loop over all tags with 1 or more yrs of missing dbh data
  tag.val <- montpines$sitetag[i]                                         # get tag no
  tag.each <- subset(montpines, montpines$sitetag==tag.val)              # Subset main data for tag i
  tag.surv <- tag.each$surv[!is.na(tag.each$surv) & tag.each$RowNum>i]   # Get non-missing years of survival after first year
  rows.wo.sz.alive$Alive[rows.wo.sz.alive$Rows==i] <- tag.surv[1]        # Store surv for 1st non-missing yr post each missed yr
  rm(tag.val,tag.each,tag.surv)
}

rows.wo.sz.alive$Alive[is.na(rows.wo.sz.alive$Alive)] <- 0  # Change NAs to 0, these are lines where missed yr was last & recorded as dead
rows.wo.sz.alive <- subset(rows.wo.sz.alive, Alive == 1)
rows.wo.sz.alive <- rows.wo.sz.alive$Rows                   # now these correspond to row numbers in main df

# Number of rows where there is missing data but tree is alive
Nindirectszcases <- length(rows.wo.sz.alive)

#rows not included in either rows.wo.sz.alive or rows.w.sz: 
included.repro.rows = c(rows.wo.sz.alive, rows.w.sz) 
rows.excluded.repro <- base::setdiff(montpines$RowNum, included.repro.rows) #setdiff finds elements in X but not in Y 
Nexcludedcases <- length(rows.excluded.repro)

## Make yr and transect numerical to use in jags as random effects ------
montpines$Year.num <- as.factor(montpines$yr) %>% as.numeric()
Year.num <- montpines$Year.num

montpines$transect.num <- as.factor(montpines$Site) %>% as.numeric()
transect.num <- montpines$transect.num

# linear index of transect-year combos
yrtranscombo=100*montpines$transect.num+montpines$Year.num

## Make dataframe w new plts (that are likely recent seedlings) for each transect & yr ------------
## only using years/trans in resurv areas 
## Make df that will hold data containing new plants

# Isolate resurvey rows of mont pines df 
resurv.montpines = montpines[montpines$resurveyarea == 1,]

# Data frame of unique transect-year combos for resurveyed areas
montpines.newPlts <- select(resurv.montpines,transect.num,Year.num) %>% unique()
colnames(montpines.newPlts) <- c("TransectNew.num", "Year.num") 

## Get first case of new plant==1 by sitetag
newPlts <- montpines %>% filter(newplt==1) %>% 
  group_by(sitetag) %>% slice(1)

#montpines[montpines$newplt == 1, ]

# Create dataframe what has the number of new plants per year-site combo
num.newPlts <- newPlts %>% # pull new plt rows from original df
  group_by("TransectNew.num"=transect.num,Year.num) %>% 
  summarise(num.newPlts=n()) # Count N of newplts per yr&trans

## Add number of new plants to df of each transect & year
montpines.newPlts <- left_join(montpines.newPlts, num.newPlts, 
                               by=c("TransectNew.num", "Year.num"))



###### These seem like relevant&misleading zeros? could be a little fishy...or maybe not? #####
## montpines.newPlts$num.newPlts[is.na(montpines.newPlts$num.newPlts)] <- 0   #Change NAs (no new plants) to zeros
montpines.newPlts$num.newPlts[montpines.newPlts$Year.num==1] <- NA         #Change new plts in year 1 to NA

ggplot(montpines.newPlts,aes(x=Year.num,y=num.newPlts,color=as.factor(TransectNew.num)))+
  geom_point()+geom_line()

## Add column so new plts in t+1 match year t. Note: variables w '1' at the end represent t+1 data
# OLD # montpines.newPlts <- montpines.newPlts %>% mutate(num.newPlts1=lead(num.newPlts)) # This needs to be done by site
transectNums = unique(montpines.newPlts$TransectNew.num)
montpines.newPlts$num.newPlts1 = NA
for(t.num in transectNums){
  
  montpines.newPlts$num.newPlts1[montpines.newPlts$TransectNew.num==t.num] = 
    lead(montpines.newPlts$num.newPlts[montpines.newPlts$TransectNew.num==t.num])
  
  }

# montpines.newPlts <- montpines.newPlts[is.na(montpines.newPlts$num.newPlts1) == FALSE,] #I am not sure why this is necessary but in erio code...
newplts <- montpines.newPlts$num.newPlts1
newplt.trans <- montpines.newPlts$TransectNew.num
newplt.yr <- montpines.newPlts$Year.num
newPltlines <- length(montpines.newPlts$TransectNew.num)

## Make a linear index of transect-year combos for new plts
newplt.yrtranscombo=100*newplt.trans+newplt.yr 

## Convert montpines new plants df into vectors. Rename columns to help keep indexing straight. 
colnames(montpines.newPlts) <- paste0('mp_newplt', colnames(montpines.newPlts))
colnames(newPlts) <- paste0('newplt_', colnames(newPlts))

mpcols <- colnames(montpines.newPlts)
for(coli in mpcols){
  tmp <- pull(montpines.newPlts,eval(as.name(coli)))
  assign(coli,tmp)
  rm(tmp)
}

mp4jags <- montpines %>% 
  select(yr,DBH,surv,
         lags,lagsrtsz,
         TempOctApr=Temp.Oct.Apr.,
         TempMaySept=Temp.May.Sept.,
         PrecipAugJuly=Precip.Aug.July,
         fog,cloud,
         TempOctApr1, TempMaySept1, 
         PrecipAugJuly1, fog1, cloud1,
         transect.num)

###########################################
# These are variables used in the prior distributions that are meant to set a climate variable to zero or not
# These values adjust the precision (inverse of variance) of the prior distribution so very large precision values here mean very small variance (i.e., forcing it to zero?)


# Big and small names are relative to uniform distibutions
big = 1e6 # This is a large precision value that will make for vary narrow priors (near zero)
small = 1e-6 # Grabbed this from vals previously used in dnorm for priors

cldmin = -small
cldmax = small

fogmin = -big
fogmax = big

tmayseptmin = -small
tmayseptmax = small

toctaprmin = -big
toctaprmax = big

precipmin = -big
precipmax = big

## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
jags.mod <- run.jags('mp_JAGSmodel_Nov2022.R', # Call to specific jags model
                     n.chains=3,
                     data=mp4jags,
                     burnin=500,
                     thin=10,
                     sample=10000,
                     adapt=100,
                     method="parallel")

# failed.jags('model')

saveRDS(jags.mod, file = 'modeloutput/zeroing_cld_sT_full.rds') 
jags.mod <- readRDS("modeloutput/zeroing_cld_sT_full.rds")


# ## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
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
