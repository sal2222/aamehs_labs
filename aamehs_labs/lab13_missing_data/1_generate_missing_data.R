# Missing Data Session
# Sebastian Rowland and Marianthi Kioumourtzoglou
# Advanced Analytic Methods for Environmental Epidemiology
# April 17, 2018
# Session 12: Missing Data


##################################################
#### Table of Contents ####
# 0: Preparation 
# 2: Create MAR Data - Confounders Only
# 3: Create MAR Data - Confounders and Exposure
# 4: Create MAR Data - Confounders, Exposure, and Outcome
#################################################

####********************
#### 0: Preparation ####
####********************

# 0a Load packages

library(readr)
library(dplyr)


# 0b Declare directories

DataPath <- "Dropbox/AAMEHS_Spring2019/Labs/Session12_MissingData/"
# 
# DataPath   <- paste0("/Users/marianthi_anna/Dropbox/AAMEHS_Spring2019/Labs/Session12_MissingData/")

####*******************
##### 1: Load Data ####
####*******************

# 1a Load data 

df <- read_csv(paste0(DataPath, "county_bmi_pm_confounders_10.csv"))

####*****************************************************
#### # 2: Create MAR missing data - Confounders Only ####
####*****************************************************
# 2: Create MAR missing data 
# we will make med.hinc, lt.HS, male.unemp have missing data 
# The probability of being missing will depend on other variables, 
# but will not depend on the value of the variable itself

# 2a Prepare dataset 
# 2a.i Remove non-continuous data 

df1 <- df %>% select(-climate.region, -AREA, -PERIMETER, -geometry)

# 2a.ii Scale the data

df1 <- scale(df1, center = TRUE)

# 2b Compute probability of missingness for each variable 

df2 <- df1 %>% 
  as.data.frame() %>%
  mutate(med.hinc.p0   = exp(-15 + -7 * female.unemp + 7 * per.black),
         lt.HS.p0      = exp(-10  + 5 * med.hval),
         male.unemp.p0 = exp(-15 + 10 * med.hval),
         avePM.idw.p0  = exp(-15 + -7 * lt.HS + 5 * med.hinc),
         aveBMI.p0     = exp(-15 + 5 * med.hinc + 2*per.white)) %>% 
  mutate(med.hinc.p = med.hinc.p0/(1+med.hinc.p0),
         lt.HS.p = lt.HS.p0/(1+lt.HS.p0),
         male.unemp.p = male.unemp.p0/(1+ male.unemp.p0),
         avePM.idw.p = avePM.idw.p0/(1+avePM.idw.p0),
         aveBMI.p = aveBMI.p0/(1+aveBMI.p0)) %>% 
  rowwise %>% 
  mutate(med.hinc.m0 = rbinom(1,1, med.hinc.p), 
         lt.HS.m0 = rbinom(1,1, lt.HS.p),
         male.unemp.m0 = rbinom(1,1, male.unemp.p),
         avePM.idw.m0 = rbinom(1,1, avePM.idw.p),
         aveBMI.m0 = rbinom(1,1, aveBMI.p)) 

# 2c Create the missingness in the variable columns   

df.m.c <- df %>% select(-AREA, -PERIMETER, -geometry, -popd)

df.m.c[df2$med.hinc.m0 == 1,]$med.hinc <- NA
df.m.c[df2$lt.HS.m0 == 1,]$lt.HS <- NA
df.m.c[df2$male.unemp.m0 == 1,]$male.unemp <- NA

# 2d Save dataframe 
df.m.c %>% write.csv(paste0(DataPath, "county_bmi_miss_confounders.csv"))


####***********************************************************
#### 3: Create MAR missing data - Confounders and Exposure ####
####***********************************************************

# 3a Create the missingness in the variable columns 

df.m.c[df2$avePM.idw.m0 == 1,]$avePM.idw <- NA   

# 3b Save dataframe 
df.m.c %>% write.csv(paste0(DataPath, "county_bmi_miss_confounders_exp.csv"))

####*******************************************************************
#### 4: Create MAR missing data - Confounders Exposure and Outcome ####
####*******************************************************************

# 4a Create the missingness in the variable columns 

df.m.c[df2$aveBMI.m0 == 1,]$aveBMI <- NA

#2d Save dataframe 
df.m.c %>% write.csv(paste0(DataPath, "county_bmi_miss_confounders_exp_out.csv"))

