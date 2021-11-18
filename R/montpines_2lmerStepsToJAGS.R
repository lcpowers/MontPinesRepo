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
library(Hmisc)
## SET WD ON LOCAL COMPUTER (WHERE JAGS SCRIPT IS LOCATED) ----------------------------------------
setwd("./R/") # should work if within project but will give error if already in this directory

## LOAD DATA --------------------------------------------------------------------------------------
montpines <- read.csv("../data/fullannual data A.csv") %>% 
  select(-1) %>% ## Get rid of first column that seems like old row numbers from elsehwere
  filter(!TagNo%in%c(0,"1(no tag)")) %>% # Get rid of what seem like extra obs. Reduces data by 3 rows
  arrange(TagNo,Site,yr)

head(montpines)

## SETTING UP NECESSARY VARIABLES -----------------------------------------------------------------
##### Setting up the jags model with lagged values

lagvals <- montpines$lags

# number of rows in data
Nallrows <- length(montpines$Site)

# number of years in data
numyears <- n_distinct(montpines$yr)

# Sites
numsites <- n_distinct(montpines$Site)

dbh <- montpines$DBH
TempOctApr <- montpines$Temp.Oct.Apr.
TempMaySept <- montpines$Temp.May.Sept.
fog <- montpines$fog
cloud <- montpines$cloud
PrecipAugJuly <- montpines$Precip.Aug.July

## Identify rows that are good dependent values (ending sizes) for surv or growth  
# I think we need a growth column?
montpines$growth <- NA
for(i in 2:nrow(montpines)) {
  
  if(montpines$TagNo[i]==montpines$TagNo[i-montpines$lags[i]]){ # if this row and lag row have the same TagNo
    montpines$growth[i] = montpines$DBH[i] - # growth = DBH of current row minus 
      montpines$DBH[i-montpines$lags[i]] # lagged DBH 
  }
}
rm(i)  

mp_summ <- montpines %>% 
  group_by(Site,yr) %>% 
  summarise(#survival=sum(surv),
            #deaths=length(surv)-sum(surv),
            n_tags = n(),
            meansz = mean(DBH,na.rm=T),
            meangrow = mean(growth,na.rm=T))

ggplot(mp_summ,aes(x=yr,y=n_tags,color=Site))+
  geom_line()+
  labs(y="Number of individuals")+
  theme_bw()

ggplot(mp_summ, aes(x=yr,y=meansz,color=Site))+
  geom_line()
  
ggplot(montpines,aes(x=DBH,y=growth,color=Site))+
  geom_point()+
  geom_smooth()+
  facet_wrap(~Site,scales='free')+
  theme_bw()

ggplot(mp_summ,aes(x=yr,y=meangrow))+
  geom_col()+
  # geom_smooth()+
  labs(y="Mean growth")+
  facet_wrap(~Site,scales='free')+
  theme_bw()

ggplot(montpines,aes(x=yr,y=growth,na.rm=T))+
  geom_point()+
  geom_smooth(method='lm')+
  facet_wrap(~Site,scales='free')+
  theme_bw()

goodrows <- which(montpines$lags>0) # This finds the rows with data
Ngoodrows <- length(goodrows) # 8625

# TO-DO: Figure out what good grow rows are for us
# goodgrowrows <- which(montpines$lags > 0)
# Ngoodgrowrows <- length(goodgrowrows) # 8625

# Vector of lag values for good rows only
goodlagvals <- montpines$lags[goodrows]

# vector of survival values
survs <- montpines$surv

##### eriogonum code has some repro code/info here. Skipping most of that for now #####
rows.w.sz <- which(!is.na(montpines$DBH)) # rows with size values
rows.wo.sz <- which(is.na(montpines$DBH)) # rows without size values
Ndirectszcases <- length(rows.w.sz)
Nindirectszcases <- length(rows.wo.sz)

## Add vector to indicate if alive or dead after missing yr(s)
montpines$RowNum <- 1:nrow(montpines) # Add a column to indicate row number
rows.wo.sz.alive <- as.data.frame(matrix(NA, # initiate data.frame
                                         nrow=length(rows.wo.sz),
                                         ncol=2))
colnames(rows.wo.sz.alive) <- c("Rows", "Alive") 
rows.wo.sz.alive$Rows <- rows.wo.sz # insert row indexes for rows without size


for (i in rows.wo.sz) {                                                  # Loop over all tags with 1 or more missing yrs                      
  tag.val <- montpines$TagNo[i] # get tag no 
  tag.each <- subset(montpines, montpines$TagNo==tag.val)                # Subset main data for tag i
  tag.surv <- tag.each$surv[!is.na(tag.each$surv) & tag.each$RowNum>i]   # Store surv for 1st non-missing yr post each missed yr
  rows.wo.sz.alive$Alive[rows.wo.sz.alive$Rows==i] <- tag.surv[1] 
  rm(tag.val,tag.each,tag.surv)
}

rows.wo.sz.alive$Alive[is.na(rows.wo.sz.alive$Alive)] <- 0  # Change NAs to 0, these are lines where missed yr was last & recorded as dead
rows.wo.sz.alive <- rows.wo.sz.alive$Alive                  # Change to vector

## At this point in eriogonum code there are lines to make transects random effects. Skipping these for now based on our conversation where Dan said it probably made more sense to treat sites as fixed effects (and no transects in the mont pines data)

## Make year a factor
# Check that year columns in mont pines data are identical
which(montpines$demoyr!=montpines$yr)

montpines$Year.num <- as.factor(montpines$demoyr) %>% as.numeric()
Year.num <- montpines$Year.num

## Set up a logical variable that is whether there is surv or grwth data in a yr (0,1): this is needed for the summing up of infl nums to predict new plts
datayesno <- rep(1,length(montpines$surv))
datayesno[which(is.na(montpines$surv)==TRUE)] <- 0 
## ------------------------------------------------------------------------------------------------

## Make dataframe w new plts (that are likely recent seedlings) for each transect & yr ------------
## Make df that will hold data containing new plants
## This code looks different than eriogonum new plants code but should to the same thing

hist(montpines$DBH[montpines$DBH<2],breaks=30)

dbh.cutoff <- 0.5

newPlts <- montpines %>% 
  group_by(TagNo) %>% # For each tag number
  slice(which.min(yr)) %>% # Find the first occurrence year
  ungroup() %>% # ungroup
  group_by(Site) %>% # For each site 
  filter(yr>min(yr) & DBH<dbh.cutoff) # filter for small plants (<0.5 dbh?) that first occurred after the site's first year

hist(newPlts$DBH,breaks=30)

num.newPlts <- newPlts %>% 
  group_by(Site,yr) %>% 
  summarise(num.newPlts=n())

## Coming back to this new plants section ##


###########################################

## RUN ASSOCIATED JAGS MODEL ----------------------------------------------------------------------
jags.mod <- run.jags('montpines_JAGSmod.R', n.chains=3, data=dats, burnin=5000, thin=10, sample=30000, adapt=500, method='parallel')

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
