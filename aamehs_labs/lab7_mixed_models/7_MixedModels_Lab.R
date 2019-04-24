# Mixed Models Lab Session
# Marianthi-Anna Kioumourtzoglou
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# March 6, 2019
# Session 7: Mixed Models 

####***********************
#### Table of Contents ####
# 0:  Preparation 
# 1:  Prepare Longitudinal Data
# 1A: Estimate ICC
# 2:  Between County Model
# 3:  Fixed Intercepts 
# 4:  Random Intercepts 
# 5:  Random Slopes and Random Intercepts 
# 6:  Compare Models
# Footnote: Cross-Sectional Model with Grouping Variable 

####********************
#### 0: Preparation ####
####********************
# 0a Install packages 

# the lme4 package provides 
# functions to create and plot 
# random intercepts and random slopes

# !! choose no, you do not want to install from sources 
#install.packages("lme4")

# The ICC package provides tools 
# to estimate the intra-class correlation 

install.packages("ICC")


# 0b Load packages

library(readr)
library(dplyr) 
library(lme4)
library(ICC) 



# 1a Readin Data 

df <- read_csv("county_bmi_pm_confounders_0509.csv")

# note that this is the similar data 
# as in the quantile regression analysis, 
# but now we have 5 years of BMI and PM2.5 data
# note that Year is a character variable 
# so in regressions Year will be a categorical variable

# 1b Convert fips to character

df <- df %>% mutate(fips = as.character(fips))

# 1c Scale and center the variables

df <- df %>% mutate(PM.scaled = scale(avePM.idw), 
                    medhinc = scale(medhinc), 
                    medhval = scale(medhval), 
                    ltHS = scale(ltHS), 
                    femaleunemp = scale(femaleunemp), 
                    maleunemp = scale(maleunemp))

####**************
#### 1A: ICC #####
####**************

# 1A.i Estimate the intra-class correlation of the outcome

BMI.ICC <- ICCest(fips, aveBMI, data = df)
BMI.ICC$ICC
BMI.ICC$LowerCI
BMI.ICC$UpperCI

# 1A.ii Estimate the intra-class correlation of the exposure

PM.ICC  <- ICCest(fips, PM.scaled, data = df)
PM.ICC$ICC
PM.ICC$LowerCI
PM.ICC$UpperCI

####******************************
#### 2: Between County Model #####
####******************************

# If we do not account for the correlation of observations within counties
# Then we must reduce our data to a single 
# observation per county in order to 
# satisfy the independence assumption 
# of linear regression models. 

# 2a Construct between - county dataset 
# Since we have multiple years, 
# we will average across years for each county

df.groups <- df %>% 
  group_by(fips) %>% 
  summarise(MeanaveBMI      = mean(aveBMI), 
            MeanPM.scaled   = mean(PM.scaled), 
            MeanMedhinc      = mean(medhinc), 
            MeanMedhval     = mean(medhval), 
            MeanltHS        = mean(ltHS), 
            MeanFemaleunemp = mean(femaleunemp), 
            MeanMaleunemp   = mean(maleunemp))

# 2b Between - county model

mod.btwn <- lm(MeanaveBMI ~ MeanPM.scaled + MeanMedhinc+ MeanMedhval+ MeanltHS + 
                 MeanFemaleunemp + MeanMaleunemp, data = df.groups)

# 2c Model summary

summary(mod.btwn)

# 2d Extract coefficients for PM 

beta.btwn.pm <- summary(mod.btwn)$coefficients[2,1]
se.btwn.pm   <- summary(mod.btwn)$coefficients[2,2]
lci.btwn.pm  <- beta.btwn.pm - 1.96*se.btwn.pm
uci.btwn.pm  <- beta.btwn.pm + 1.96*se.btwn.pm

coef.btwn.pm <- paste0(round(beta.btwn.pm,3), " (95% CI: ", round(lci.btwn.pm,3), ", ", round(uci.btwn.pm,3), ")")
coef.btwn.pm 

# 2e Translate coefficients to 1 unit increase in PM2.5
# When we scaled the data, 
# we divided each variable by its standard deviation 
# We can translate the coefficients to the original units 
# by multiplying by the standard deviation of the data 

# 2e.i Compute standard deviation

sd.pm <- sd(df$avePM.idw)

# 2e.ii Multiply the estimate and 95% CI 

beta.btwn.pm1 <- sd.pm * beta.btwn.pm
lci.btwn.pm1 <- sd.pm * lci.btwn.pm
uci.btwn.pm1 <- sd.pm * uci.btwn.pm

coef.btwn.pm1 <- paste0(round(beta.btwn.pm1,3), " (95% CI: ", round(lci.btwn.pm1,3), ", ", round(uci.btwn.pm1,3), ")")
coef.btwn.pm1 

####****************************
#### 3:  Fixed Intercepts  #####
####****************************

# 3a Model with fixed intercept 
# Within subject effect estimate
# The model will estimate 235 intercepts  
# one for each county. 

mod.fix.int <- lm(aveBMI ~ PM.scaled +   # exposure
                    Year +               # confounders
                    fips,                # fixed intercept for fips
                  data = df)

# 3b Coefficients of model 

coef.table.fix.int <- summary(mod.fix.int)$coefficient
dim(coef.table.fix.int) 

# the model estimated coefficients for 
# the overall intercept, 1 exposure,
# 4 confounders (one per year) (except year =05 is reference category)
# and 234 intercepts for counties 
# one county was treated as the baseline county
# as with any categorical variable

# Just for fun:

coef.table.fix.int

# 3c Estimates for PM2.5 

beta.FI.pm <- coef.table.fix.int[2,1] 
se.FI.pm   <- coef.table.fix.int[2,2]
lci.FI.pm  <- beta.FI.pm - 1.96*se.FI.pm
uci.FI.pm  <- beta.FI.pm + 1.96*se.FI.pm

coef.FI.pm <- paste0(round(beta.FI.pm,3), " (95% CI: ", round(lci.FI.pm,3), ", ", round(uci.FI.pm,3), ")")
coef.FI.pm 

# 3d Compare with between-county model

coef.btwn.pm 

# why now different than above?
# PM.scaled    0.14944    0.03457   4.323 2.22e-05 ***
# what assumption do each of these models have?
# (think ICC)

####***************************
#### 4: Random Intercept  #####
####***************************

# 4a Model with Random Intercept 
# the (1|fips) signifies a random intercept for each county

mod.ran.int <- lmer(aveBMI ~ PM.scaled +              # exposure
                      Year + medhinc+ medhval+ ltHS + # confounders
                      femaleunemp + maleunemp +       # confounders
                      (1|fips),                       # random intercept for each county
                    data = df) 

# 4b Model Summary 

summary(mod.ran.int)
summary(mod.ran.int)$coefficients

beta.RI.pm <- summary(mod.ran.int)$coefficients[2,1]
se.RI.pm   <- summary(mod.ran.int)$coefficients[2,2]
lci.RI.pm  <- beta.RI.pm - 1.96*se.RI.pm
uci.RI.pm  <- beta.RI.pm + 1.96*se.RI.pm


coef.RI.pm <- paste0(round(beta.RI.pm,3), " (95% CI: ", round(lci.RI.pm,3), ", ", round(uci.RI.pm,3), ")")
coef.RI.pm 

# 4c Compare to previous models 

coef.FI.pm
# and
coef.btwn.pm

####**********************************************
#### 5: Random Slopes and Random Intercepts  #####
####**********************************************

# 5a Model with Random Intercept and Random slope
# the (avePM.idw|fips) indicates a random slope of pm-bmi 
# for each county. 
# the (1|fips) signifies a random intercept for each county


mod.ran.int.slope <- lmer(aveBMI ~ PM.scaled +              # exposure
                            Year + medhinc+ medhval+ ltHS + # confounders
                            femaleunemp + maleunemp +       # confounders
                            (1 + PM.scaled|fips),           # random terms 
                          data = df)

# 5b Model Summary 

summary(mod.ran.int.slope)

beta.RE.pm <- summary(mod.ran.int.slope)$coefficients[2,1]
se.RE.pm   <- summary(mod.ran.int.slope)$coefficients[2,2]
lci.RE.pm  <- beta.RE.pm - 1.96*se.RE.pm
uci.RE.pm  <- beta.RE.pm + 1.96*se.RE.pm

coef.RE.pm <- paste0(round(beta.RE.pm,3), " (95% CI: ", round(lci.RE.pm,3), ", ", round(uci.RE.pm,3), ")")
coef.RE.pm 

# 5c Plot Random Slopes 

ran.eff <- ranef(mod.ran.int.slope)$fips

# 5c.ii Density plot 

plot(density(ran.eff$PM.scaled))

# 5c.iii Density of total PM association in each county. 

tot.pm <- ran.eff$PM.scaled + beta.RE.pm
plot(density(tot.pm))
lines(density(ran.eff$PM.scaled), col = "red")

####*************************
#### 6: Compare Models  #####
####*************************
# 6a Compare Model with fixed intercept and model with random intercept
# Likelihood Ratio Test

anova( mod.ran.int, mod.fix.int)

# 6b Compare model with random intercept and model with random intercept and random slope 
# Likelihood Ratio Test

anova( mod.ran.int, mod.ran.int.slope)


####***************************************************************
#### Footnote: Cross-sectional Model with Grouping Variable   #####
####***************************************************************
# A Readin Data 

df.10 <- read_csv(paste0(DataPath, "county_bmi_pm_confounders_10.csv"))

# B Create State Data 
df.10 <- df.10 %>% mutate(state = substr(fips, 0,2))

# C Cross-sectional Model 
# random intercept for state 

mod.cs <- lmer(aveBMI ~ avePM.idw + per.black + per.latinx + per.asnam+ 
                 medhinc + medhval +ltHS + femaleunemp + 
                 (1|state),
               data = df.10)

# D Model summary 

summary(mod.cs)

# E Model Coefficients 

coef.table <- summary(mod.cs)$coefficients

