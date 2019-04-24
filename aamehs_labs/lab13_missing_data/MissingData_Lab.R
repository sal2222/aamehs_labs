# Missing Data Session
# Sebastian Rowland and Marianthi Kioumourtzoglou
# Advanced Analytic Methods for Environmental Epidemiology
# April 17, 2018
# Session 12: Missing Data

##################################################
#### Table of Contents ####
# 0: Preparation 

## A: Only Confounder Data Missing
# A1: Load MAR Data - Confounders Only
# A2: Visualize Missingness
# A3: Impute Missing Data 
# A4: Fit Model on Imputed Data
# A5: Pool Models
# A6: Compare to Model from Original Data

## B: Confounder and Exposure Data Missing 
# B1: Load MAR Data - Confounders and Exposure
# B2: Visualize Missingness
# B3: Impute Missing Data 
# B4: Fit Model on Imputed Data
# B5: Pool Models
# B6: Compare to Model from Original Data

## C: Confounder, Exposure, and Outcome Data Missing 
# C1: Load MAR Data - Confounders, Exposure, and Outcome
# C2: Visualize Missingness
# C3: Impute Missing Data 
# C4: Fit Model on Imputed Data
# C5: Pool Models
# C6: Compare to Model from Original Data

####********************
#### 0: Preparation ####
####********************
# 0a Install packages 

# install.packages("mice")
# install.packages("VIM")

# 0b Load packages

library(readr)
library(dplyr)
library(mice)
library(VIM)

# 0c Declare directories

DataPath <- "Dropbox/AAMEHS_Spring2019/Labs/Session12_MissingData/"
# 
# DataPath   <- paste0("C:/Users/mk3961/Dropbox/AAMEHS_Spring2019/Labs/Session12_MissingData/")
# 
# 

# 0d Turn off scientific notation

options(scipen=999)

###############################################################

##*************************************##
#### A: Only Confounder Data Missing ####
##*************************************##

####********************
##### A1: Load Data ####
####********************

# 1a Load data 

df.m.c <- read_csv(paste0(DataPath, "county_bmi_miss_confounders.csv"))

# 1b Review data

head(df.m.c)

# 1c Remove empty column 

df.m.c <- df.m.c %>% select(-X1)

####******************************
#### A2: Visualize Missingness ####
####******************************
# 2a Count number of missing for each variable

summary(df.m.c)

sum(is.na(df.m.c$med.hinc))
sum(is.na(df.m.c$lt.HS))
sum(is.na(df.m.c$male.unemp))
sum(is.na(df.m.c$avePM.idw))
sum(is.na(df.m.c$aveBMI))

# 2b Review pattern of missingness

# We will isolate just the variables with missing data so it is easier to see them 
df.miss <- df.m.c %>% select(med.hinc, lt.HS, male.unemp)

# 2b.i Visualize with mice package

md.pattern(df.miss)

# 2b.ii Visualize with VIM package 

aggr_plot <- aggr(df.miss, 
                  col=c('navyblue','red'), 
                  numbers=TRUE,
                  sortVars=TRUE, 
                  labels=names(df.miss), 
                  cex.axis=.7, 
                  gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))

# 2c Visual missingness of two variables 

df.miss.var <- df.miss %>% dplyr::select(med.hinc, lt.HS)
dev.off()
marginplot(df.miss.var, xlab = "Median Household Income", ylab = "Percent Less than HS")

####*****************************
#### A3: Impute Missing Data ####
####*****************************

# 3a Create imputed datasets

df.mi.c <- mice(df.m.c,     # dataset
              m = 5,        # number of imputed datasets to create
              maxit = 50,   # number of iterations to create prediction models 
              meth = 'pmm', # method for imputation- depends on the class of missing variables
                            # pmm usually performs best for continuous variables
              seed = 500)   # our seed

# 3b Review imputed datasets

summary(df.mi.c)

# 3c Review imputed values for a single variable 

df.mi.c$imp$med.hinc

# 3d Compare density plot of imputed data versus original data 

densityplot(df.mi.c)

####***********************************
#### A4: Fit Model on Imputed Data ####
####***********************************

# 4A Fit models

mod1.c <- with(df.mi.c, lm(aveBMI ~ avePM.idw + med.hinc + per.black + per.latinx + per.asnam + 
                         med.hval + lt.HS + female.unemp + male.unemp + climate.region))

####*********************
#### A5: Pool Models ####
####*********************

# 5a Model summary 

pool(mod1.c)

# ubar - average within-imputation variance (sample variance)
# b - between-imputation variance

# riv - relative increase in variance due to missingness
# lambda - proportion of variance due to missingness 
# fmi- fraction of missing information

# we can directly compute the standard error of Beta 
# from the two sources of variance
# Var = se^2
# se = sqrt(ubar = ((m+ 1)/m)* b))
sqrt(0.0007955+ ((5 + 1)/5) * 0.000000335992)

# 5b Model summary II 
# summary() will compute the standard error for us

mod.pooled.c <- summary(pool(mod1.c))
mod.pooled.c

####*********************************************
#### A6: Compare to Model from Original Data ####
####*********************************************

# 6a Readin original, complete data 

df <- read_csv(paste0(DataPath, "county_bmi_pm_confounders_10.csv"))

# 6b Create model with original data

mod.df <- summary(lm(aveBMI ~ avePM.idw + med.hinc + per.black + per.latinx + per.asnam + 
     med.hval + lt.HS + female.unemp + male.unemp + climate.region, data = df))

# 6c Compare coefficient estimates

mod.df$coefficients[2:6,1]
mod.pooled.c[2:6,1]

# 6c Compare standard errors 

mod.df$coefficients[2:6,2]
mod.pooled.c[2:6,2]


#
#
##**********************************************##
#### B: Confounder and Exposure Data Missing  ####
##**********************************************##

####********************
##### B1: Load Data ####
####********************

# 1a Load data 

df.m.ce <- read_csv(paste0(DataPath, "county_bmi_miss_confounders_exp.csv"))

# 1b Remove empty column 

df.m.ce <- df.m.ce %>% select(-X1)

####******************************
#### B2: Visualize Missingness ####
####******************************
# 2a Count number of missing for each variable

# summary(df.m.ce)

sum(is.na(df.m.ce$med.hinc))
sum(is.na(df.m.ce$lt.HS))
sum(is.na(df.m.ce$male.unemp))
sum(is.na(df.m.ce$avePM.idw))
sum(is.na(df.m.ce$aveBMI))


# 2b Review pattern of missingness

# We will isolate just the variables with missing data so it is easier to see them 
df.miss <- df.m.ce %>% select(med.hinc, lt.HS, male.unemp, avePM.idw)

# 2b.i Visualize with VIM package 

aggr_plot <- aggr(df.miss, 
                  col=c('navyblue','red'), 
                  numbers=TRUE,
                  sortVars=TRUE, 
                  labels=names(df.miss), 
                  cex.axis=.7, 
                  gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))

# 2c Visual missingness of two variables 

# marginplot(select(df.miss, avePM.idw, lt.HS))

####*****************************
#### B3: Impute Missing Data ####
####*****************************

# 3a Create imputed datasets

df.mi.ce <- mice(df.m.ce,      # dataset
              m = 5,       # number of imputed datasets to create
              maxit = 50,   # number of iterations to create prediction models 
              meth = 'pmm', # method for imputation- depends on the class of missing variables
                            # pmm usually performs best for continuous variables
              seed = 500,   # our seed
              printFlag = FALSE) # silence iteration names
              
# 3b Review imputed datasets

# summary(df.mi.ce)

# 3c Review imputed values for a single variable 

# df.mi.ce$imp$avePM.idw

# 3d Compare density plot of imputed data versus original data 

densityplot(df.mi.ce)

# 3e Plot convergence trace plots
# Plots show the mean and std of the inputed values 
# where each color is a imputed dataset (m = 10)
# 
# plot(df.mi.ce, c("male.unemp", "lt.HS"))
# 
# plot(df.mi.ce, c("lt.HS", "avePM.idw"))

####***********************************
#### B4: Fit Model on Imputed Data ####
####***********************************

# 4A Fit models

mod1.ce <- with(df.mi.ce, lm(aveBMI ~ avePM.idw + med.hinc + per.black + per.latinx + per.asnam + 
                         med.hval + lt.HS + female.unemp + male.unemp + climate.region))

####*********************
#### B5: Pool Models ####
####*********************

# 5a Model summary 

mod.pooled.ce <- summary(pool(mod1.ce))
mod.pooled.ce

####*********************************************
#### B6: Compare to Model from Original Data ####
####*********************************************

# 6a Create model with original data

mod.df <- summary(lm(aveBMI ~ avePM.idw + med.hinc + per.black + per.latinx + per.asnam + 
                       med.hval + lt.HS + female.unemp + male.unemp + climate.region, data = df))

# 6b Compare coefficient estimates

mod.df$coefficients[2:6,1]
mod.pooled.c[2:6,1]
mod.pooled.ce[2:6,1]

# 6c Compare standard errors 

mod.df$coefficients[2:6,2]
mod.pooled.c[2:6,2]
mod.pooled.ce[2:6,2]



##*******************************************************##
#### C: Confounder, Exposure and Outcome Data Missing  ####
##*******************************************************##

####********************
##### C1: Load Data ####
####********************

# 1a Load data 

df.m.ceo <- read_csv(paste0(DataPath, "county_bmi_miss_confounders_exp_out.csv"))

# 1b Review data

# head(df.m.ceo)

# 1c Remove empty column 

df.m.ceo <- df.m.ceo %>% select(-X1)

####*******************************
#### C2: Visualize Missingness ####
####*******************************
# 2a Count number of missing for each variable

# summary(df.m.ceo)

sum(is.na(df.m.ceo$med.hinc))
sum(is.na(df.m.ceo$lt.HS))
sum(is.na(df.m.ceo$male.unemp))
sum(is.na(df.m.ceo$avePM.idw))
sum(is.na(df.m.ceo$aveBMI))


# 2b Review pattern of missingness

# We will isolate just the variables with missing data so it is easier to see them 
df.miss <- df.m.ceo %>% select(med.hinc, lt.HS, male.unemp, avePM.idw, aveBMI)

# 2b.i Visualize with mice package

# md.pattern(df.miss)

# 2b.ii Visualize with VIM package 

aggr_plot <- aggr(df.miss, 
                  col=c('navyblue','red'), 
                  numbers=TRUE,
                  sortVars=TRUE, 
                  labels=names(df.miss), 
                  cex.axis=.7, 
                  gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))

# 2c Visual missingness of two variables 

# marginplot(select(df.miss, med.hinc, aveBMI))

####*****************************
#### C3: Impute Missing Data ####
####*****************************

# 3a Create imputed datasets

df.mi.ceo <- mice(df.m.ceo,          # dataset
              m = 5,       # number of imputed datasets to create
              maxit = 50,   # number of iterations to create prediction models 
              meth = 'pmm', # method for imputation- depends on the class of missing variables
                            # pmm usually performs best for continuous variables
              seed = 500,   # our seed
              printFlag = FALSE) # silence iteration names
# 3b Review imputed datasets

# summary(df.mi.ceo)

# 3c Review imputed values for a single variable 

# df.mi.ceo$imp$aveBMI

# 3d Compare density plot of imputed data versus original data 

densityplot(df.mi.ceo)

####***********************************
#### C4: Fit Model on Imputed Data ####
####***********************************

# 4A Fit models

mod1.ceo <- with(df.mi.ceo, lm(aveBMI ~ avePM.idw + med.hinc + per.black + per.latinx + per.asnam + 
                         med.hval + lt.HS + female.unemp + male.unemp + climate.region))

####*********************
#### C5: Pool Models ####
####*********************

# 5a Model summary 

mod.pooled.ceo <- summary(pool(mod1.ceo))
mod.pooled.ceo

####*********************************************
#### C6: Compare to Model from Original Data ####
####*********************************************

# 6a Create model with original data

mod.df <- summary(lm(aveBMI ~ avePM.idw + med.hinc + per.black + per.latinx + per.asnam + 
                       med.hval + lt.HS + female.unemp + male.unemp + climate.region, data = df))

# 6b Compare coefficient estimates

mod.df$coefficients[2:6,1]
mod.pooled.c[2:6,1]
mod.pooled.ce[2:6,1]
mod.pooled.ceo[2:6,1]

# 6c Compare standard errors 

mod.df$coefficients[2:6,2]
mod.pooled.c[2:6,2]
mod.pooled.ce[2:6,2]
mod.pooled.ceo[2:6,2]







# Bonus:  Plot convergence trace plots
# Plots show the mean and std of the inputed values 
# where each color is a imputed dataset (m = 10)

plot(df.mi.c, c("male.unemp", "lt.HS", "med.hinc"))

