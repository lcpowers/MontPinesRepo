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
setwd("./R/") # should work if within project but will give error if already in 'R' folder

## LOAD DATA --------------------------------------------------------------------------------------
montpines <- read.csv("../data/fullannual data A.csv") %>% 
  select(-1) %>% ## Get rid of first column that seems like old row numbers from elsehwere
  filter(!TagNo%in%c(0,"1(no tag)")) %>% # Get rid of what seem like extra obs. Reduces data by 3 rows
  arrange(TagNo,Site,yr)

# head(montpines)

montpines$growth <- NA # Size this year minus size last measured year - skips lag>1 years
montpines$propgrowth <- NA # (size this year minus size last year)/size last measured year skips lag >1 years
montpines$allgrowth <- NA # Growth for all years -- deals with lag >1 years
montpines$allpropgrowth <- NA # Proportional growth for all years -- deals with lag >1 years

for(i in 2:nrow(montpines)) {
  
  if(montpines$TagNo[i]==montpines$TagNo[i-montpines$lags[i]]){ # if this row and lag row have the same TagNo
    montpines$growth[i] = montpines$DBH[i] - montpines$DBH[i-montpines$lags[i]] # absolute growth = DBH of current row - lagged DBH
    montpines$propgrowth[i] = montpines$growth[i]/montpines$DBH[i-montpines$lags[i]]
  }
  if(montpines$lags[i]==1) montpines$allgrowth[i] = montpines$growth[i] # If there is data for two consecutive years (lag=1), don't divide growth
  
  if(montpines$lags[i]>1){ # If the lag is greater than 1,evenely distribute growth between all skipped years
  l = montpines$lags[i]
  growthvec = rep(montpines$growth[i]/l,l)
  montpines$allgrowth[(i-l+1):i] = growthvec
  propgrowthvec = rep(montpines$proprgrowth[i]/l,l)
  montpines$allproprgrowth[(i-l+1):i] = propgrowthvec
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
lagvals <- montpines$lags #Full list of lags or of -1 for first observation rows
goodrows <- which(montpines$lags>0) # This finds the rows with data
goodgrowrows <- which(montpines$lagsrtsz > 0) # This finds rows with good growth data

goodlagvals <- montpines$lags[goodrows]
Ncases <- length(goodrows)
Ngrowcases <- length(goodgrowrows)
Survs <- montpines$surv
dbh <- montpines$DBH


##### eriogonum code has some repro code/info here. Skipping most of that for now #####
rows.w.sz <- which(!is.na(montpines$DBH)) # rows with size values
rows.wo.sz <- which(is.na(montpines$DBH)) # rows without size values
Ndirectszcases <- length(rows.w.sz)
Nindirectszcases <- length(rows.wo.sz)

## create vector to indicate if alive or dead after missing yr(s)
montpines$RowNum <- 1:nrow(montpines) # Add a column to indicate row number
rows.wo.sz.alive <- as.data.frame(matrix(NA, # initiate data.frame
                                         nrow=length(rows.wo.sz),
                                         ncol=2))
colnames(rows.wo.sz.alive) <- c("Rows", "Alive") 
rows.wo.sz.alive$Rows <- rows.wo.sz # insert row indexes for rows without size

for (i in rows.wo.sz) {                                                  # Loop over all tags with 1 or more yrs of missing dbh data               
  tag.val <- montpines$TagNo[i]                                          # get tag no 
  tag.each <- subset(montpines, montpines$TagNo==tag.val)                # Subset main data for tag i
  tag.surv <- tag.each$surv[!is.na(tag.each$surv) & tag.each$RowNum>i]   # Get non-missing years of survival after first year
  rows.wo.sz.alive$Alive[rows.wo.sz.alive$Rows==i] <- tag.surv[1]        # Store surv for 1st non-missing yr post each missed yr
  rm(tag.val,tag.each,tag.surv)
}

rows.wo.sz.alive$Alive[is.na(rows.wo.sz.alive$Alive)] <- 0  # Change NAs to 0, these are lines where missed yr was last & recorded as dead
rows.wo.sz.alive <- rows.wo.sz.alive$Alive                  # Change to vector

## At this point in eriogonum code there are lines to make transects random effects. Skipping these for now based on our conversation where Dan said it probably made more sense to treat sites as fixed effects (and no transects in the mont pines data)

## Make year numerical to use in jags as random effects 
# Check that year columns in mont pines data are identical
montpines$Year.num <- as.factor(montpines$demoyr) %>% as.numeric()
Year.num <- montpines$Year.num

##### Should we do this with site? #####

## Set up a logical variable that is whether there is surv or grwth data in a yr (0,1): this is needed for the summing up of infl nums to predict new plts
datayesno <- rep(1,length(montpines$surv))
datayesno[which(is.na(montpines$surv)==TRUE)] <- 0 
## ------------------------------------------------------------------------------------------------

## Make dataframe w new plts (that are likely recent seedlings) for each transect & yr ------------
## Make df that will hold data containing new plants
## This code looks different than eriogonum new plants code but should to the same thing
years <- unique(montpines$Year.num)
years <- years[order(years)]

# hist(montpines$DBH[montpines$DBH<2],breaks=30)
dbh.cutoff <- 0.5
newPlts <- montpines %>% 
  group_by(TagNo) %>% # For each tag number
  slice(which.min(yr)) %>% # Find the first occurrence year for each tag and take the first obs
  ungroup() %>% # ungroup
  group_by(Site) %>% # For each site 
  filter(yr>min(yr) & DBH<dbh.cutoff) # filter for small plants (<0.5 dbh?) that first occurred after the site's first year
# hist(newPlts$DBH,breaks=30)

num.newPlts <- newPlts %>% 
  group_by(Site,yr) %>% 
  summarise(num.newPlts=n())

##### TO-DO: Come back to this new plants section #####

mp4jags <- montpines %>% 
  select(yr,DBH,surv,
         lags,lagsrtsz,
         TempOctApr=Temp.Oct.Apr.,
         TempMaySept=Temp.May.Sept.,
         PrecipAugJuly=Precip.Aug.July,
         fog,cloud,lagsrtsz)

###########################################

## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
jags.mod <- run.jags('montpines_JAGSmodel.R', n.chains=1, data=mp4jags, burnin=5000, thin=10, sample=30000, adapt=500)

failed.jags('model')

#save(jags.mod, file='erbr_JAGSmod_c3t10s20b5_210406.rdata')
# saveRDS(jags.mod, "erbr_JAGSmod_c3t10s30b5_210509.rds")
## ------------------------------------------------------------------------------------------------

# 
# ## LOOK AT MODEL OUTPUT ---------------------------------------------------------------------------
# jags.mod <- readRDS("erbr_JAGSmod_c3t10s30b10_noYRE_210613.rds")
# summary(jags.mod)
# plot(jags.mod)
# summ.mod <- summary(jags.mod)
# tail(summ.mod[,1:5], n=32)
# gelman.diag(jags.mod, confidence = 0.95, transform=FALSE)
# 
# chains <- jags.mod$mcmc
# chains <- bind_rows(lapply(chains, as.data.frame))
# colMeds <- apply(chains,2,median)
# colSDs <- apply(chains,2,sd)
# 
# 
# ## Look at correlation b/w params 
# chains.1 <- chains %>% dplyr::select(!contains(c("randomeffect", "precision")))
# chains.1 <- chains.1 %>% dplyr::select(!c(deviance, resid.sum.sq))
# 
# cor.chains <- cor(chains.1)
# corrplot(cor.chains, method="circle", type="lower")


## ** Make bar graph comparing median param ests & 80% CIs b/w diff datasets **
## ------------------------------------------------------------------------------------------------  


