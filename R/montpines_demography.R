## This is a rough script for turning output of JAGS model into 
# matrix model 

montpines <- read.csv("../data/fullannual data C.csv") %>% 
  select(-1)  ## Get rid of first column that seems like old row numbers from elsewhere
montpines <- subset(montpines, montpines$demoyr < 2015) #wonky years

#JAGS model output 
output <- readRDS("jags_mod_output_Apr18_22.rds")
chains <- output$mcmc
chains <- bind_rows(lapply(chains, as.data.frame))
colMeds <- as.data.frame(apply(chains,2,median))
Medians <- t(cbind.data.frame(colMeds$`apply(chains, 2, median)`))
colnames(Medians) = rownames(colMeds)
Medians = as.data.frame(Medians)
attach(Medians)

#get range of weather years, mu = avg, coldest, and hottest: 
mufog = mean(montpines$fog)
muTempMaySept = mean(montpines$Temp.May.Sept.)
muPrecipAugJuly = mean(montpines$Precip.Aug.July)

coldyrfog = max(montpines$fog)
coldyrTempMaySept = min(montpines$TempMaySept)
coldyrPrecipAugJuly = max(montpines$Precip.Aug.July)

hotyrfog = min(montpines$fog)
hotyrTempMaySept = max(montpines$Temp.May.Sept.)
hotyrPrecipAugJuly = min(montpines$Precip.Aug.July)

# first thing is to get estimates of survival, grwth for a range of DBHs 
DBHs = seq(min(montpines$DBH, na.rm = T), max(montpines$DBH, na.rm = T), 5)

# Survival estimates, avg years. This is just using same form of the model 
# in JAGS right now and plugging in estimates of each coefficient/value for each DBH
# done for each transect in each type of weather yr (5 transects * 3 weather types)

surv_estimates_trans1_mu = 1/(1+exp(-(surv_intercept+
                              DBHs*surv_dbhCoef+
                              surv_TempMaySeptCoef*muTempMaySept+
                              surv_PrecipAugJulyCoef*muPrecipAugJuly+
                              surv_fogCoef*mufog+
                              `surv_Transect_randomeffect[1]`)))

surv_estimates_trans2_mu = 1/(1+exp(-(surv_intercept+
                                        DBHs*surv_dbhCoef+
                                        surv_TempMaySeptCoef*muTempMaySept+
                                        surv_PrecipAugJulyCoef*muPrecipAugJuly+
                                        surv_fogCoef*mufog+
                                        `surv_Transect_randomeffect[2]`)))

surv_estimates_trans3_mu = 1/(1+exp(-(surv_intercept+
                                        DBHs*surv_dbhCoef+
                                        surv_TempMaySeptCoef*muTempMaySept+
                                        surv_PrecipAugJulyCoef*muPrecipAugJuly+
                                        surv_fogCoef*mufog+
                                        `surv_Transect_randomeffect[3]`)))

surv_estimates_trans4_mu = 1/(1+exp(-(surv_intercept+
                                        DBHs*surv_dbhCoef+
                                        surv_TempMaySeptCoef*muTempMaySept+
                                        surv_PrecipAugJulyCoef*muPrecipAugJuly+
                                        surv_fogCoef*mufog+
                                        `surv_Transect_randomeffect[4]`)))

surv_estimates_trans5_mu = 1/(1+exp(-(surv_intercept+
                                        DBHs*surv_dbhCoef+
                                        surv_TempMaySeptCoef*muTempMaySept+
                                        surv_PrecipAugJulyCoef*muPrecipAugJuly+
                                        surv_fogCoef*mufog+
                                        `surv_Transect_randomeffect[5]`)))

#Surv estimates, hot dry years: 
surv_estimates_trans1_hotdry = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                        surv_TempMaySeptCoef*hotyrTempMaySept+
                                        surv_PrecipAugJulyCoef*hotyrPrecipAugJuly+
                                        surv_fogCoef*hotyrfog+`surv_Transect_randomeffect[1]`)))

surv_estimates_trans2_hotdry = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                            surv_TempMaySeptCoef*hotyrTempMaySept+
                                            surv_PrecipAugJulyCoef*hotyrPrecipAugJuly+
                                            surv_fogCoef*hotyrfog+`surv_Transect_randomeffect[2]`)))

surv_estimates_trans3_hotdry = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                            surv_TempMaySeptCoef*hotyrTempMaySept+
                                            surv_PrecipAugJulyCoef*hotyrPrecipAugJuly+
                                            surv_fogCoef*hotyrfog+`surv_Transect_randomeffect[3]`)))

surv_estimates_trans4_hotdry = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                            surv_TempMaySeptCoef*hotyrTempMaySept+
                                            surv_PrecipAugJulyCoef*hotyrPrecipAugJuly+
                                            surv_fogCoef*hotyrfog+`surv_Transect_randomeffect[4]`)))

surv_estimates_trans5_hotdry = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                            surv_TempMaySeptCoef*hotyrTempMaySept+
                                            surv_PrecipAugJulyCoef*hotyrPrecipAugJuly+
                                            surv_fogCoef*hotyrfog+`surv_Transect_randomeffect[5]`)))

#cold wet yr 
surv_estimates_trans1_cold = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                          surv_TempMaySeptCoef*coldyrTempMaySept+
                                          surv_PrecipAugJulyCoef*coldyrPrecipAugJuly+
                                          surv_fogCoef*coldyrfog+`surv_Transect_randomeffect[1]`)))

surv_estimates_trans2_cold = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                          surv_TempMaySeptCoef*coldyrTempMaySept+
                                          surv_PrecipAugJulyCoef*coldyrPrecipAugJuly+
                                          surv_fogCoef*coldyrfog+`surv_Transect_randomeffect[2]`)))

surv_estimates_trans3_cold = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                          surv_TempMaySeptCoef*coldyrTempMaySept+
                                          surv_PrecipAugJulyCoef*coldyrPrecipAugJuly+
                                          surv_fogCoef*coldyrfog+`surv_Transect_randomeffect[3]`)))

surv_estimates_trans4_cold = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                          surv_TempMaySeptCoef*coldyrTempMaySept+
                                          surv_PrecipAugJulyCoef*coldyrPrecipAugJuly+
                                          surv_fogCoef*coldyrfog+`surv_Transect_randomeffect[4]`)))

surv_estimates_trans5_cold = 1/(1+exp(-(surv_intercept+DBHs*surv_dbhCoef+
                                            surv_TempMaySeptCoef*coldyrTempMaySept+
                                            surv_PrecipAugJulyCoef*coldyrPrecipAugJuly+
                                            surv_fogCoef*coldyrfog+`surv_Transect_randomeffect[5]`)))


# sorry this is so ugly. but here we go, just turning these estimates into a df I can easily use 
surv_outs = cbind.data.frame(DBHs, surv_estimates_trans1_mu, surv_estimates_trans2_mu, surv_estimates_trans3_mu, 
                             surv_estimates_trans4_mu, surv_estimates_trans5_mu, 
                             surv_estimates_trans1_cold, surv_estimates_trans2_cold, surv_estimates_trans3_cold, 
                             surv_estimates_trans4_cold, surv_estimates_trans5_cold, 
                             surv_estimates_trans1_hotdry, surv_estimates_trans2_hotdry, surv_estimates_trans3_hotdry, 
                             surv_estimates_trans4_hotdry, surv_estimates_trans5_hotdry)

colnames(surv_outs) = c('DBHs', 't1mu.surv','t2mu.surv','t3mu.surv','t4mu.surv','t5mu.surv',
                        't1cold.surv','t2cold.surv','t3cold.surv','t4cold.surv','t5cold.surv',
                        't1hot.surv','t2hot.surv','t3hot.surv','t4hot.surv','t5hot.surv')

###### Repro estimates ##### this is potentially really not right.  
repro_sizes <- reproDBHcoef*DBHs + repro_intercept 
repro_sizes[1] =  0.01 #let's keep it positive  

#repro is a poisson process in JAGS, so we'll use our lambda values 
repros = qpois(.5, lambda = log(repro_sizes)) 
repros[1] = 0
repros= cbind.data.frame(DBHs, repros)

## Growth estimates ### 
# same as before, using the model in same format as JAGS and solving for estimated size at t+1 
# for each DBH, given each transect*climate effect 
# We need to come back and add in correct grwthvar_dbhCoef which wasn't monitored on last model run.  
grwthvar_dbhCoef = 0.07245125
sz_estimates_trans1_mu <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                 grwth_TempMaySeptCoef*muTempMaySept + 
                                grwth_PrecipAugJulyCoef*muPrecipAugJuly +
                                  grwth_fogCoef*mufog +                 
                                  `grwth_Transect_randomeffect[1]`, 
                                exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans2_mu <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                         grwth_TempMaySeptCoef*muTempMaySept + 
                                                         grwth_PrecipAugJulyCoef*muPrecipAugJuly +
                                                         grwth_fogCoef*mufog +                 
                                                         `grwth_Transect_randomeffect[2]` , 
                                                         exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans3_mu <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                         grwth_TempMaySeptCoef*muTempMaySept + 
                                                         grwth_PrecipAugJulyCoef*muPrecipAugJuly +
                                                         grwth_fogCoef*mufog +                 
                                                         `grwth_Transect_randomeffect[3]`, 
                                                         exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans4_mu <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                         grwth_TempMaySeptCoef*muTempMaySept + 
                                                         grwth_PrecipAugJulyCoef*muPrecipAugJuly +
                                                         grwth_fogCoef*mufog +                 
                                                         `grwth_Transect_randomeffect[4]`, 
                                                         exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans5_mu <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                         grwth_TempMaySeptCoef*muTempMaySept + 
                                                         grwth_PrecipAugJulyCoef*muPrecipAugJuly +
                                                         grwth_fogCoef*mufog +                 
                                                         `grwth_Transect_randomeffect[5]`, 
                                                         exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}

sz_estimates_trans1_hot <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                              grwth_TempMaySeptCoef*hotyrTempMaySept + 
                                                              grwth_PrecipAugJulyCoef*hotyrPrecipAugJuly +
                                                              grwth_fogCoef*hotyrfog +                 
                                                              `grwth_Transect_randomeffect[1]`, 
                                                            exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans2_hot <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                               grwth_TempMaySeptCoef*hotyrTempMaySept + 
                                                               grwth_PrecipAugJulyCoef*hotyrPrecipAugJuly +
                                                               grwth_fogCoef*hotyrfog +                 
                                                               `grwth_Transect_randomeffect[2]`, 
                                                             exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans3_hot <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                               grwth_TempMaySeptCoef*hotyrTempMaySept + 
                                                               grwth_PrecipAugJulyCoef*hotyrPrecipAugJuly +
                                                               grwth_fogCoef*hotyrfog +                 
                                                               `grwth_Transect_randomeffect[3]`, 
                                                             exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans4_hot <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                               grwth_TempMaySeptCoef*hotyrTempMaySept + 
                                                               grwth_PrecipAugJulyCoef*hotyrPrecipAugJuly +
                                                               grwth_fogCoef*hotyrfog +                 
                                                               `grwth_Transect_randomeffect[4]`, 
                                                             exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans5_hot <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                               grwth_TempMaySeptCoef*hotyrTempMaySept + 
                                                               grwth_PrecipAugJulyCoef*hotyrPrecipAugJuly +
                                                               grwth_fogCoef*hotyrfog +                 
                                                               `grwth_Transect_randomeffect[5]`, 
                                                             exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}

sz_estimates_trans1_cold <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                               grwth_TempMaySeptCoef*coldyrTempMaySept + 
                                                               grwth_PrecipAugJulyCoef*coldyrPrecipAugJuly +
                                                               grwth_fogCoef*coldyrfog +                 
                                                               `grwth_Transect_randomeffect[1]`, 
                                                             exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans2_cold <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                                grwth_TempMaySeptCoef*coldyrTempMaySept + 
                                                                grwth_PrecipAugJulyCoef*coldyrPrecipAugJuly +
                                                                grwth_fogCoef*coldyrfog +                 
                                                                `grwth_Transect_randomeffect[2]`, 
                                                              exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans3_cold <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                                grwth_TempMaySeptCoef*coldyrTempMaySept + 
                                                                grwth_PrecipAugJulyCoef*coldyrPrecipAugJuly +
                                                                grwth_fogCoef*coldyrfog +                 
                                                                `grwth_Transect_randomeffect[3]`, 
                                                              exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans4_cold <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                                grwth_TempMaySeptCoef*coldyrTempMaySept + 
                                                                grwth_PrecipAugJulyCoef*coldyrPrecipAugJuly +
                                                                grwth_fogCoef*coldyrfog +                 
                                                                `grwth_Transect_randomeffect[4]`, 
                                                              exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}
sz_estimates_trans5_cold <- function(currentDBH) {return(list(grwth_intercept + grwth_dbhCoef*currentDBH + 
                                                                grwth_TempMaySeptCoef*coldyrTempMaySept + 
                                                                grwth_PrecipAugJulyCoef*coldyrPrecipAugJuly +
                                                                grwth_fogCoef*coldyrfog +                 
                                                                `grwth_Transect_randomeffect[5]`, 
                                                              exp(grwthvar_intercept + grwthvar_dbhCoef*currentDBH)))}

#trying to copy Dan to turn growth into cdf fns for each bin:
nbin = length(DBHs) #number to integrate over? or # of bins? 
minsize = min(montpines$DBH, na.rm =T) 
maxsize = (max(montpines$DBH, na.rm = T))
vecbin = seq(minsize + .01, maxsize +.01, 1)

prob_growth_outs = matrix(nrow = length(DBHs), ncol = 15)

#put all growth functions in a list 
grwth_functions = c(sz_estimates_trans1_mu,sz_estimates_trans2_mu,
                    sz_estimates_trans3_mu,sz_estimates_trans4_mu,sz_estimates_trans5_mu,
                    sz_estimates_trans1_cold,sz_estimates_trans2_cold,
                    sz_estimates_trans3_cold,sz_estimates_trans4_cold,sz_estimates_trans5_cold,
                    sz_estimates_trans1_hot,sz_estimates_trans2_hot,
                    sz_estimates_trans3_hot,sz_estimates_trans4_hot,sz_estimates_trans5_hot)

## loop to calculate sz transition probabilities for all 15 scenarios 

for (f in (1:length(grwth_functions))){
  
  fn = grwth_functions[[f]]
  
for (ss in 1:(nbin-1)) {
  
  sizes_in_bin = seq(DBHs[ss]+.01, DBHs[ss+1], .05)
  
  #for each size, what is prob of moving to next bin? 
  sizes_t1 = qnorm(.5, mean = fn(sizes_in_bin)[[1]], 
        sd = sqrt(fn(sizes_in_bin)[[2]]))
  
  #if those are new growing sizes, how many do we think we expect to 
  # 'escape' into the next bin? 
  next_bin = DBHs[ss+1]
  prop_next_bin = length(which(sizes_t1 >= next_bin))/length(sizes_t1)
  
  #save prop of trees growing into next bin
  prob_growth_outs[ss,f] = prop_next_bin
}
}

grwth_outs = cbind.data.frame(DBHs, prob_growth_outs)

colnames(grwth_outs) = c('DBHs', 't1mu','t2mu','t3mu','t4mu','t5mu',
                         't1cold','t2cold','t3cold','t4cold','t5cold',
                         't1hot','t2hot','t3hot','t4hot','t5hot')

grwth_outs[21,2:16] = 1 #make sure we have 100% prob of staying in last size class

########## Time for matrix construction!!! ############## 
grwth_outs # these is our growth rates 
repros #repro values (no random effect of transect or climate in JAGS model yet)
surv_outs #survival rates 


## skeleton of matrices: 

#        0-5              5-10                 10-15 .........   100+
#0-5   Repro*Surv +      Repro * Surv           Repro * Surv     Repro * Surv
#     prob(don't grow)

#5-10 Surv*prob_grow    Surv*prob(dontgrow)       0                      0

#10-15   0              Surv*prob(grow)          Surv*prob(dontgrow)

#....    0                0                      Surv * prob(grow)
 
#100+    0                0

matrices = list(NULL)  

years = c('mu','cold','hot') #this order matters - lines up with surv_outs
transects= seq(1,5,1)
i = 0

for (yy in years){
    
    for (t in transects){
    
    i = i+1
    mx_name = paste(yy,t,sep = "_") #is this correct? 
    print(mx_name)
    mx = matrix(0,ncol= length(DBHs), nrow = length(DBHs))   #initialize mx
    mx[1,] = repros$repros * surv_outs[,i+1]                   #this fills in top row w/ Surv * Repro
    
    #fill in diagonals with prob of NOT growing * surv
    diags = surv_outs[,i+1]*(1-grwth_outs[,i+1])
    diag(mx) = diags
    
    #fill in that next diag with prob of growing * surv
    offdiag = surv_outs[,i+1]*grwth_outs[,i+1]    
    for (cols in 1:20) {mx[cols + 1, cols] = offdiag[cols]}
    #mx[1,1] = mx[1,1] + (1-grwth_outs[1,i+1])*surv_outs[1,i+1]            #for first stage need to adjust by prob of not growing
    
    #save final mx
    matrices[[i]] =  mx
  }
}

### lambda variable setup and plotting ### 

lam_fun = function(x) {return((eigen(x)$values[1]))}
lambdas = lapply(matrices, FUN = lam_fun)
lambdas = unlist(lambdas)
lambdas = as.data.frame(Re(lambdas))

lambdas$transects = rep(seq(1,5,1),3)
lambdas$yrs = rep(years, each = 5)
lambdas$site = rep(unique(montpines$Site), 3)
colnames(lambdas) = c('lambda', 'transect','year_type', 'site')
lambdas$year_type = factor(lambdas$year_type, levels = c('cold','hot','mu'))

#make a plot just of lambdas
ggplot(lambdas) + geom_point(aes(site, lambda, col = year_type), size = 6) + 
  geom_hline(yintercept = 1, lty = 2) + theme_bw(base_size = 14) + 
  scale_color_manual(values=c("#56B4E9","#FF0000","#999999")) + 
  ylim(.5, 1.5)


#make a plot just to show adult tree biomass? 
summed_DBHs = aggregate(montpines$DBH, 
                        by = list(montpines$Site, montpines$demoyr), 
                        FUN = sum, na.rm = T)

colnames(summed_DBHs) = c('Site','yrnum','summedDBH')

# lets lump the 2009-2010 MSC data together, because seems like they did half survey in each year. 
# 2010 would've been "correct" survey year. 
summed_DBHs[(summed_DBHs$Site == "MSC" & summed_DBHs$yrnum == "2010"),]$summedDBH <- 2061
summed_DBHs[(summed_DBHs$Site == "MSC" & summed_DBHs$yrnum == "2009"),]$summedDBH <- NA

mp = left_join(montpines, summed_DBHs, by = c('Site', 'demoyr' = 'yrnum'))

mpyrs = subset(mp, mp$summedDBH > 80) #more sleuthing - some odd super low values still lingering?

ggplot(mpyrs) + 
  geom_point(aes(demoyr, summedDBH, col = Site), size = 5) + theme_bw() + 
  facet_wrap(vars(Site), scales = "free_y") + 
  labs(col = 'transect') + xlab("demoyr") + ylab("summedDBH") + 
  geom_line(aes(demoyr, summedDBH, col = Site)) + 
  scale_color_viridis(discrete = TRUE)
  

# show the climate data along w numbers thru time
ssc = subset(mp, mp$Site == "SSC")
ggplot(ssc) + 
  geom_point(aes(demoyr, summedDBH)) + theme_bw(base_size = 20) +
  xlab("Year") + ylab("summed DBH") + 
  geom_line(aes(demoyr, summedDBH)) + 
  geom_point(aes(demoyr, cloud/10), col = 'blue') + 
  geom_line(aes(demoyr, cloud/10), col = 'blue')+ 
  scale_y_continuous(name = 'summedDBH', 
                     sec.axis = sec_axis(~./10, name ='cloud'))

ggplot(mp) + 
  geom_point(aes(demoyr, summedDBH)) + theme_bw() + 
  facet_wrap(vars(Site), scales = "free_y") + 
  xlab("demoyr") + ylab("summedDBH") + 
  geom_line(aes(demoyr, summedDBH)) + 
  geom_point(aes(demoyr, cloud/10), col = 'blue') + 
  geom_line(aes(demoyr, cloud/10), col = 'blue')+ 
  scale_y_continuous(name = 'summedDBH', 
                     sec.axis = sec_axis(~./10, name ='cloud'))

#add secondary axis
ggplot(mpyrs) + 
  geom_point(aes(demoyr, summedDBH)) + theme_bw() + 
  facet_wrap(vars(Site), scales = "free_y") + 
  xlab("demoyr") + ylab("summedDBH") + 
  geom_line(aes(demoyr, summedDBH)) + 
  geom_line(aes(demoyr,Temp.Oct.Apr.*100), col = 'blue') + 
  geom_point(aes(demoyr, Temp.Oct.Apr.*100), col = 'blue') + 
  scale_y_continuous(name = 'summedDBH', 
                     sec.axis = sec_axis(~.*10, name ='Winter Temp'))

#add secondary axis
ggplot(mpyrs) + 
  geom_point(aes(demoyr, summedDBH)) + theme_bw() + 
  facet_wrap(vars(Site), scales = "free_y") + 
  xlab("demoyr") + ylab("summedDBH") + 
  geom_line(aes(demoyr, summedDBH)) + 
  geom_line(aes(demoyr,Precip.Aug.July), col = 'blue') + 
  geom_point(aes(demoyr, Precip.Aug.July), col = 'blue') + 
  scale_y_continuous(name = 'summedDBH', 
                     sec.axis = sec_axis(~., name ='Summer Precip'))+ 
  labs(legend = c('DBHs','summer precip'))
  


## add a simulation - say next 100 years are mix of all years 
## now lets say they're mostly hot, some average, few cold..? 
## we're going to start with our # of trees that we have now 

montpines$szbins = round_any(montpines$DBH, 5, f = floor)

#get counts of current Ns in each size bin in 2014 @ each transect
Ns = aggregate(montpines$szbins, by = list(montpines$Site), FUN = count)
lastyr = subset(montpines, montpines$demoyr == 2014) 
counts = lastyr[,c("Site","szbins")] 
counts$tally = rep(1, length(counts$Site))
counts = subset(counts, counts$szbins != "NA")
starting_distribution = aggregate(counts$tally, by = list(counts$Site,counts$szbins), FUN = sum)
colnames(starting_distribution) = c('Site', 'SzClass','N')
all_bins = rep(seq(0, 100, 5),5)
all_sites = rep(unique(montpines$Site), 21)
zeros = cbind.data.frame(all_bins, rep(0, length(all_bins)))
zeros = cbind.data.frame(zeros, all_sites)
Nis = merge(starting_distribution, zeros, by.x = c('SzClass','Site'), by.y = c('all_bins','all_sites'), all = T) 
Nis[is.na(Nis)] <- 0

#reorder & start simulation 
Nis = arrange(Nis, Nis$Site, Nis$SzClass)
Ni1 = as.vector(Nis$N[1:21])
Ni2 = as.vector(Nis$N[1:21])
years = 10
sims= 100
totalNs = matrix(NA, nrow = 10, ncol = 500) 

colnum = 0
  
 for (t in transects){
   
  for (s in (1:sims)){
     
     #set N0 and restart
      colnum = colnum + 1
      N0 = as.vector(Nis$N[t:(t+20)]) 
      Nt = N0
  
    for (yy in (1:years)){
    
    #grab correct indices from matrices and select 1 of the 3 yr types
    mx_is = c(t, t+5, t+10)
    climate = sample(mx_is, 1)
    
    Nt1 = matrices[[climate]] %*% Nt
    Nt = Nt1 #new Nt for next time 
    
    # save count 
    totalNs[yy,colnum] = round(sum(Nt), 0)
    
    }
   }
 } 


#organize output from simulation 
totalNs = as.data.frame(totalNs)
totalNs$t = seq(1,years)
m = reshape2::melt(totalNs, id.vars = 't')
m$SITE = rep(unique(montpines$Site), each = 1000)
m$variable = NULL 

#take median in each year
medians = aggregate(m$value, by = list(m$SITE, m$t), FUN = median)
colnames(medians) = c('Site', 't', 'median')
ggplot(medians) + geom_line(aes(t, median, col = Site)) + 
  theme_bw(base_size = 14) + 
  scale_color_viridis(discrete = T) +
  xlab('years') + ylab('N')
