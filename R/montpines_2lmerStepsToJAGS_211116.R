## Doak Lab Monterey Pines Project
## Data preparation and call to JAGS script
## Modify data and assign variables needed for the JAGS model with data lags (missing years)
## Info about associated JAGS script later
##
#### This file is the equivalent of <erbr_2lmerStepsToJAGS_210504.R>
#### I've going through that file line-by-line, bringing over the code that seems relevant to us. 
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
library(Hmisc)
## SET WD ON LOCAL COMPUTER (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------
setwd("./R/")

## LOAD DATA --------------------------------------------------------------------------------------
montpines <- read.csv("../fullannual data A.csv") %>% 
  select(-1) %>% ## Get rid of first column that seems like old row numbers from elsehwere
  filter(!TagNo%in%c(0,"1(no tag)")) %>% 
  arrange(TagNo,Site,yr)

head(montpines)

## SETTING UP NECESSARY VARIABLES -----------------------------------------------------------------
##### Setting up the jags model with lagged values
Nallrows <- length(montpines$Site)

# Years in data
numyears <- n_distinct(montpines$yr)

# Sites
numsites <- n_distinct(montpines$Site)
  
## Identify rows that are good dependent values (ending sizes) for surv or growth  
# I think we need a growth column...

montpines <- montpines %>% 
  group_by(TagNo) %>% 
  mutate(growth = lag(x = DBH,n=max(lags,1)))

montpines$growth <- NA

for(i in 2:nrow(montpines)) {
  
  if(montpines$TagNo[i]==montpines$TagNo[i-montpines$lags[i]]){ # if this row and lag row are the same individual
    montpines$growth[i] = montpines$DBH[i] - # DBH of current row minus
      montpines$DBH[i-montpines$lags[i]] # lagged DBH 
  }
  
}
rm(i)  
  
goodrows <- which(montpines$lags>0) # This finds the rows with data
Ngoodrows <- length(goodrows)

# Vector of lag values for good rows only
goodlagvals <- montpines$lags[goodrows]

# vector of survival values
survs <- montpines$surv


##### April's code has repro info here. Skipping for now #####
rows.w.sz <- which(!is.na(montpines$DBH)) # rows with size values
rows.wo.sz <- which(is.na(montpines$DBH)) # rows without size values


## Add vector to indicate if alive or dead after missing yr(s)

montpines$RowNum <- 1:nrow(montpines) # Add a column to indicate row number

rows.wo.sz.alive <- as.data.frame(matrix(NA, nrow=length(rows.wo.sz), ncol=2))
colnames(rows.wo.sz.alive) <- c("Rows", "Alive")
rows.wo.sz.alive$Rows <- rows.wo.sz


hist(montpines$DBH[montpines$DBH<3],breaks=30)
# Fitler

newPlts <- montpines %>% 
  group_by(TagNo) %>% # For each tag number
  slice(which.min(yr)) %>% # Find the first occurrence year
  ungroup() %>% # ungroup
  group_by(Site) %>% # For each site, 
  filter(yr>min(yr) & DBH<0.5) # filter for small plants that occurred after the first year





  
