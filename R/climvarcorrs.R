## Looking at correlations between climate variables
## Generally, summer and winter T, and fog/cloud are strongly correlated


library(tidyverse)
library(corrplot)
library(ggfortify)
setwd("./R/") 

# This reads in raw data and gets unique climate vars
# clim <- read.csv("../data/fullannual data C.csv") %>% 
#   select(site=Site,yr=demoyr,wT=Temp.Oct.Apr.,sT=Temp.May.Sept.,fog,cloud,ppt=Precip.Aug.July) %>% 
#   unique() %>% 
#   na.omit()

# Need to run data prep scipt through growth loop (~line 117) for this to work
clim = montpines %>% 
  select(yr,site=Site,wT=Temp.Oct.Apr.,sT=Temp.May.Sept.,fog,cloud,ppt=Precip.Aug.July,growth,anngrowth=allgrowth) %>% 
  unique() %>% 
  na.omit()

ggplot(clim,aes(x=yr,y=cloud,color=site))+
  geom_line()+
  geom_point()

corrplot(cor(clim[,3:7]),type="lower",diag = T)

climpca = prcomp(clim[,3:7],scale=T)

autoplot(climpca,data=clim,loadings=T,loadings.label=T)+
  theme_bw()

clim2 = select(clim,-c(site,contains("growth"))) %>% unique() %>% na.omit()
corrplot(cor(clim2[,2:6]),type="lower",diag = T)

clim2pca = prcomp(clim2[,2:6],scale=T)



autoplot(clim2pca,data=clim2,loadings=T,loadings.label=T)+
  theme_bw()
