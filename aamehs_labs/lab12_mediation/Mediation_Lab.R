# Mediation Session
# Sebastian Rowland and Marianthi Kioumourtzoglou
# Advanced Analytic Methods for Environmental Epidemiology
# April 17, 2018
# Session 11: Mediation
# Thanks to Linda Valeri for help developing this lab

##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Load Data 
# 2: Traditional Mediation Analysis: Difference Method
# 3: Traditional Mediation Analysis: Product Method
# 4: Mediation Package - no Interaction
# 5: Assess Interaction
# 6: Mediation Analysis With Interaction


####********************
#### 0: Preparation ####
####********************
# 0a Install packages 

install.packages("mediation")

# 0b Load packages

library(readr)
library(dplyr)
library(mediation)


# 
###############################################################
####*******************
##### 1: Load Data ####
####*******************
# 1a Load data 

df <- read_csv("county_bmi_pm_confounders_10.csv")

# 1b Review data

head(df)

summary(df)

# 1c Rescale income variable 
# income is measured in dollars, but we are more interested in changes per thousand dollars 

df <- df %>% mutate(med.hinc = med.hinc/1000)

####***********************************************************
##### 2: Traditional Mediation Analysis: Difference Method ####
####***********************************************************

# 2a  Estimate total effect

# 2a.i Fit total effects model 
# here the coefficient for median household income 
# captures the total effect between median household income 
# and averageBMI, including pathways through PM2.5
# note that we do not inclue the potential mediator (avePM.idw)
# in this model 

TE.fit  <- lm(aveBMI ~ med.hinc + per.black + per.latinx + per.asnam + 
                 med.hval + lt.HS + female.unemp + male.unemp + climate.region, 
              data = df, na.action = na.omit)

# 2a.ii Extract total effect

summary(TE.fit)

TE <- summary(TE.fit)$coeff[2,1]
TE

# 2b Estimate direct effect

# 2b.i Fit direct effects model w/o interaction
# here the coefficient for median household income 
# captures the direct effect between median household income  
# and average BMI, through pathways that do not include increased PM2.5 exposure
# note that we include the potential mediator in the model 
# not that we do not allow for interaction between median household income and PM2.5

DE.fit.noint <- lm(aveBMI ~ med.hinc + avePM.idw + per.black + per.latinx + per.asnam + 
                      med.hval + lt.HS + female.unemp + male.unemp + climate.region, 
                   data = df, na.action = na.omit)

# 2a.i.i Extract direct effect
summary(DE.fit.noint)

DE <- summary(DE.fit.noint)$coeff[2,1]
DE

# 2c Compute indirect effect via difference method 
# we subtract the direct effect from the total difference

IE.diff <- TE - DE 
IE.diff

####********************************************************
##### 3: Traditional Mediation Analysis: Product Method ####
####********************************************************

# 3a Estimate effect of exposure on mediator

# 3a.i Fit the model for the mediator 
# here the coefficient for med.hinc 
# captures the effect between the county's median household income
# and the county's annual PM2.5 

MED.fit  <- lm(avePM.idw ~  med.hinc + per.black + per.latinx + per.asnam + 
                 med.hval + lt.HS + female.unemp + male.unemp + climate.region, 
               data = df, na.action = na.omit)

# 3a.ii Extract effect of median household income on annual PM

summary(MED.fit)

exp_to_med <- summary(MED.fit)$coeff[2,1]
exp_to_med 

# 3b Estimate effect of mediator on outcome 

# 3b.i Fit direct effects model w/o interaction
# Now we are interested in estimating the effects between the mediator 
# and the outcome 
# not that we do not allow for interaction between median household income and PM2.5
# we assume that the impact of PM2.5 on average BMI is consistent across counties. 

DE.fit.noint <- lm(aveBMI ~ med.hinc + avePM.idw + per.black + per.latinx + per.asnam + 
                     med.hval + lt.HS + female.unemp + male.unemp + climate.region, 
                   data = df, na.action = na.omit)

# 3a.ii Extract effect of mediator on outcome 
summary(DE.fit.noint)

med_to_out <- summary(DE.fit.noint)$coeff[3,1]
med_to_out

# 3c Compute indirect effect via product method 
# we multiply the effect of the exposure on the potential mediator 
# with the effect of the mediator on the outcome

IE.prod <- exp_to_med * med_to_out 

# 3d Compare methods 
# the methods yield the same results

IE.diff
IE.prod


####**********************************************
##### 4: Mediation Package - no Interaction   ####
####**********************************************

# 4a Run mediation analysis under assumption of no interation 

set.seed(7) ## set seed because Imai et al. approach is based on sampling

mediation.noint <- mediate(MED.fit,        # exposure -> mediator model 
                           DE.fit.noint,   # direct effects model
                           treat = "med.hinc", 
                           mediator = "avePM.idw", 
                           control.value = 40,  
                           treat.value = 55,    
                           sims = 1000)

# 4b Model Summary 

summary(mediation.noint)

## NIE = 0.01283
## NDE = -0.03529
## TE  = -0.0224

## this is the same as both the product and difference methods 
## but scaled by the difference between the control and treament levels 

mediation.noint$d1 /15
mediation.noint$z1 /15
mediation.noint$tau.coef /15

# 4c Extract effect estimates and CI 

mediation.results.df <- data.frame(
  Estimate_Name = c("Indirect Effect",         "Direct Effect",             "Total Effect",                   "Proportion Mediated"),
  Estimate      = c(mediation.noint$d1,         mediation.noint$z1,         mediation.noint$tau.coef,         mediation.noint$n1),
  Lower_95_CI   = c(mediation.noint$d1.ci[[1]], mediation.noint$z1.ci[[1]], mediation.noint$tau.ci[[1]], mediation.noint$n1.ci[[1]]),
  Upper_95_CI   = c(mediation.noint$d1.ci[[2]], mediation.noint$z1.ci[[2]], mediation.noint$tau.ci[[2]], mediation.noint$n1.ci[[2]]),
  PValue        = c(mediation.noint$d1.p,       mediation.noint$z1.p,       mediation.noint$tau.p,       mediation.noint$n1.p))

mediation.results.df

####*****************************
##### 5: Assess Interaction  ####
####*****************************

# 5a Fit model with interaction

DE.fit.int <- lm(aveBMI ~ med.hinc * avePM.idw + 
                   per.black +  per.latinx + per.asnam + 
                   med.hval + lt.HS + female.unemp + male.unemp + climate.region, 
                 data = df, na.action = na.omit)

# 5b Model summary

summary(DE.fit.int)

# we see very weak suggestion of interaction 
# the effect of annual PM is lower among counties with high median household income

# 5c 95% Confidence intervals for interaction 

coeff.est <- summary(DE.fit.int)$coeff[19,1]
se <- summary(DE.fit.int)$coeff[19,2]
interaction.term <- c(coeff.est, coeff.est - 1.96* se, coeff.est + 1.96* se )
interaction.term



####**********************************************
##### 6: Mediation Analysis With Interaction  ####
####**********************************************

# 6a Fit mediation model

set.seed(10) ## set seed because Imai et al. approach is based on sampling

mediation.int <- mediate(MED.fit, 
                         DE.fit.int,         # Here we are using the outcome model with interaction
                         treat = "med.hinc", 
                         mediator = "avePM.idw", 
                         control.value = 40,  
                         treat.value = 55,     
                         sims = 1000)

#6b Model Summary 

summary(mediation.int)

## TNIE = ACME (treated) = 0.01179
## PNDE = ADE  (control) = -0.02069
## TE  = -0.00890

# 6c Plot indirect, direct, and total effects

plot(mediation.int)
