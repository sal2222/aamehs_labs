# DLNM Lab Session
# Marianthi-Anna Kioumourtzoglou
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# Feburary 27, 2019
# Session 6: DLNM 

####***********************
#### Table of Contents ####
# 0: Preparation 
# 1: Prepare NY Daily Data
# 2: Individual Models of Each Lag
# 3: All Lags Mutually Adjusted
# 4: Linear DLM 
# 5: Nonlinear DLNM
# FootNote: Examining Different Lag DF for Linear Model

####********************
#### 0: Preparation ####
####********************
# 0a Install packages 

# the dlnm package provides a set of functions to 
# create and plot 
# distributed lag nonlinear models 

#install.packages("dlnm")

# 0b Load packages

library(readr)
library(dplyr) 
library(lubridate)
library(splines)
library(ggplot2)
library(dlnm)


####*******************************
#### 1: Prepare NY Daily Data #####
####*******************************
# 1a Readin data 

df <- read_csv("nyc_time_series.csv")

# 1b Look at data structure 
# Data Dictionary 
# Each row is a unique day
# dateD is the date. 
# DailyPM - population-weighted average of AQS estimates across NYC 
# DailyTemp - population-weighted average of temperature in Celcius from NLDAS 
# Daily RH - population-weighted average of relative humidity in percentage from NLDAS 
# DailyCVDin - number of in-patient hospitalizations for cardiovascular disease

# 1c Convert date column to datetime format 

df <- df %>% mutate(Date = parse_date_time(dateD, "mdy")) %>% 
             select(-dateD)

# 1d Extract the day of the week

df <- df %>% mutate(DayofWeek = as.character(wday(Date, label = TRUE)))
head(df$DayofWeek)

# 1e Arrange data by date 

df <- df %>% arrange(Date)

# 1f Compute the number of years 

yr_num <- length(unique(year(df$Date)))

####**************************************
#### 2: Individual Models of Each Lag ####
####**************************************

# 2a Assign lagged PM, temperature, RH
# we will use the lag() function to assign the lagged exposures 
# but this is not the only possible way to assign lags

# pm_0 is the same variable as DailyPM 
# it is not necesary to create a new variable in order to do dlnm 
# but it helps illustrate the lags. 

df <- df %>% mutate(
  pm_0 = dailyPM,
  pm_1 = lag(dailyPM, 1), 
  pm_2 = lag(dailyPM, 2), 
  pm_3 = lag(dailyPM, 3), 
  pm_4 = lag(dailyPM, 4), 
  pm_5 = lag(dailyPM, 5), 
  pm_6 = lag(dailyPM, 6))

# repeat for temp and RH 

df <- df %>% mutate(
  temp_0 = dailyTemp,
  temp_1 = lag(dailyTemp, 1), 
  temp_2 = lag(dailyTemp, 2), 
  temp_3 = lag(dailyTemp, 3), 
  temp_4 = lag(dailyTemp, 4), 
  temp_5 = lag(dailyTemp, 5), 
  temp_6 = lag(dailyTemp, 6), 
  rh_0 = dailyRH,
  rh_1 = lag(dailyRH, 1), 
  rh_2 = lag(dailyRH, 2), 
  rh_3 = lag(dailyRH, 3), 
  rh_4 = lag(dailyRH, 4), 
  rh_5 = lag(dailyRH, 5), 
  rh_6 = lag(dailyRH, 6))

# note that the first 6 days of the study will have 
# at least one lag missing 

# 2b Look at data structure 

head(df)
View(df)

# 2c Construct the Time Series models 
# one for each lag 
# we will assume a linear dose-response relationship 
# between PM and CVD admission rate 
# and a nonlinear relationship for temperature and relative humidity 
# indicator variables for day of the week
# again, we will use a spline term for secular trends 
# with 4 degrees of freedom per year

mod.lag0 <- glm(dailyCVDin ~               # outcome
                  pm_0 +                   # exposure
                  ns(temp_0, df = 4) +     # nonlinear term for Temp
                  ns(rh_0, df = 3) +       # nonlinear term for RH
                  DayofWeek +              # categorical variable for DoW
                  ns(Date, df = 4*yr_num), # nonlinear term for secular trend
                family = "quasipoisson",   # distribution family
                data = df)

mod.lag1 <- glm(dailyCVDin ~ pm_1 +  ns(temp_1, df = 4) + 
                  ns(rh_1, df = 3) + DayofWeek + ns(Date, df = 4*yr_num), 
                family = "quasipoisson", data = df)

mod.lag2 <- glm(dailyCVDin ~ pm_2 +  ns(temp_2, df = 4) + 
                  ns(rh_2, df = 3) + DayofWeek + ns(Date, df = 4*yr_num), 
                family = "quasipoisson", data = df)

mod.lag3 <- glm(dailyCVDin ~ pm_3 +  ns(temp_3, df = 4) + 
                  ns(rh_3, df = 3) + DayofWeek + ns(Date, df = 4*yr_num), 
                family = "quasipoisson", data = df)

mod.lag4 <- glm(dailyCVDin ~ pm_4 +  ns(temp_4, df = 4) + 
                  ns(rh_4, df = 3) + DayofWeek + ns(Date, df = 4*yr_num), 
                family = "quasipoisson", data = df)

mod.lag5 <- glm(dailyCVDin ~ pm_5 +  ns(temp_5, df = 4) + 
                  ns(rh_5, df = 3) + DayofWeek + ns(Date, df = 4*yr_num), 
                family = "quasipoisson", data = df)

mod.lag6 <- glm(dailyCVDin ~ pm_6 +  ns(temp_6, df = 4) + 
                  ns(rh_6, df = 3) + DayofWeek + ns(Date, df = 4*yr_num), 
                family = "quasipoisson", data = df)

# 2d Plot the models' results 

# 2d.i Extract the coefficients from each model 
# In this code I have condensed several steps: 
# we are summarizing each model, then extracting
# the coefficient and se of the second row (pm_0)
# of the coefficients matrix

coeff.lag0 <- summary(mod.lag0)$coefficients[2, 1:2]
coeff.lag1 <- summary(mod.lag1)$coefficients[2, 1:2]
coeff.lag2 <- summary(mod.lag2)$coefficients[2, 1:2]
coeff.lag3 <- summary(mod.lag3)$coefficients[2, 1:2]
coeff.lag4 <- summary(mod.lag4)$coefficients[2, 1:2]
coeff.lag5 <- summary(mod.lag5)$coefficients[2, 1:2]
coeff.lag6 <- summary(mod.lag6)$coefficients[2, 1:2]

# 2d.ii Create dataframe of model results

coeff.table <- rbind(coeff.lag0, coeff.lag1, coeff.lag2, coeff.lag3,
                     coeff.lag4, coeff.lag5, coeff.lag6)

coeff.table <- as.data.frame(coeff.table, stringsAsFactors = FALSE)

# 2d.iii Set names for dataframe

names(coeff.table) <- c("coeff", "se")

# 2d.iv Compute confidence intervals 

coeff.table <- coeff.table %>% 
  mutate(lci = coeff - 1.96 * se, 
         uci = coeff + 1.96 * se)

# 2d.v Exponeniate terms 
# since the link function is log() 
# when we exponantiate, we translate our coefficients 
# from log(RateRatio)'s
# to Rate Ratios 
# we will also multiply by 10 to compute the 
# the rate ratio per 10 ug/m3 increase in daily PM

coeff.table <- coeff.table %>% 
  mutate(coeff.exp = exp(10*coeff), 
         lci.exp = exp(10*lci), 
         uci.exp = exp(10*uci))

# 2d.vi Create model names 

coeff.table        <- coeff.table %>% 
  mutate(ModelName = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 5", "Lag 6"))

# 2d.vii Plot
# note that this is the same code 
# that we used in the quantile regression session
# to create the forest plots

fp.individual.lags  <- ggplot(
                              # defines what dataset ggplot will use
                              data = coeff.table,     
                              # aes() defines which variables the geoms will use   
                              aes( # defines variable for the x axis
                                x = ModelName,  
                                # defines the variable for the point along the y axis
                                y = coeff.exp,      
                                # defines the lower bound of the confidence interval
                                ymin = lci.exp,     
                                # define the upper bound of the confidence interval 
                                ymax = uci.exp)) +  
                              # creates a point (y) with line defined by ymin and ymax
                              geom_pointrange() +   
                              # creates lines with bars, i.e. here the CIs
                              geom_errorbar() +      
                              # add a dashed line at y=0
                              geom_hline(aes(yintercept = 1.0), lty = 2) +
                              # labels for axes
                              xlab("Model Name") +    
                              ylab(expression("RR per 10 "*mu*"g/m"^3*" PM"[2.5]~" (95% CI)"))

fp.individual.lags

# Class Question ### 
### Do we have any issues in our interpretation of the coefficients? # 

####***********************************
#### 3: All Lags Mutually Adjusted ####
####***********************************
# We will assume each lag has a linear dose-response relationship
# In this model we will adjust for 7-day mean temperature 
# and 7-day mean relative humidity 

# 3a Compute 7-day mean temp and RH

df$meanTemp.7day <- rowMeans(df[,15:21])
df$meanRH.7day   <- rowMeans(df[,22:28])

# 3b Construct the model 

mod.lag0_6 <- glm(dailyCVDin ~ pm_0 + pm_1 + pm_2 +
                    pm_3 + pm_4 + pm_5 + pm_6+
                    ns(meanTemp.7day, df = 4) + 
                    ns(meanRH.7day, df = 3) + DayofWeek + ns(Date, df = 4*yr_num), 
                  family = "quasipoisson", data = df)

# 3c Model Summary 
summary(mod.lag0_6)

# 3d Plot model 

# 3d.i Extract Coefficients 

coeff.table <- summary(mod.lag0_6)$coefficients[2:8, 1:2]

coeff.table <- as.data.frame(coeff.table)

# 3d.ii Set names for dataframe

names(coeff.table) <- c("coeff", "se")

# 3d.iii Compute confidence intervals 

coeff.table <- coeff.table %>% 
  mutate(lci = coeff - 1.96 * se, 
         uci = coeff + 1.96 * se)

# 3d.iv Exponeniate terms 
# since the link function is log() 
# when we exponantiate, we translate our coefficients 
# from log(RateRatio)'s
# to Rate Ratios 
# we will also multiply by 10 to compute the 
# the rate ratio per 10 ug/m3 increase in daily PM

coeff.table <- coeff.table %>% 
  mutate(coeff.exp = exp(10*coeff), 
         lci.exp = exp(10*lci), 
         uci.exp = exp(10*uci))

# 3d.v Create Lag names 

coeff.table <- coeff.table %>% 
  mutate(Lag = c("Lag 0", "Lag 1", "Lag 2", "Lag 3", "Lag 4", "Lag 5", "Lag 6"))

# 3d.vi Plot

fp.mutual.adjusted.lags <- ggplot(data=coeff.table,     # defines what dataset we are using
                                  aes(x=Lag,            # defines variable for the x axis
                                      y=coeff.exp,      # defines the variable for the point along the y axis
                                      ymin=lci.exp,     # defines the lower bound of the confidence interval
                                      ymax=uci.exp)) +  # define the upper bound of the confidence interval   
                                  geom_pointrange() +   # creates a point (y) with line defined by ymin and ymax        
                                  geom_errorbar()+      # creates lines with bars
                                  geom_hline(aes(yintercept=1), lty=2) + # add a dashed line at y=1 
                                  xlab("Model Name") +                   # labels for axes
                                  ylab(expression("RR per 10 "*mu*"g/m"^3*" PM"[2.5]~" (95% CI)"))

fp.mutual.adjusted.lags

####*******************
#### 4: Linear DLM ####
####*******************

# 4a Autocorrelation plot
# acf() computes the autocorrelation 
# for each lag

acf(df$dailyPM, lag.max = 7)

# the acf appears to increase at the end, and that seems to be due to 
# a day of the week pattern in PM 

acf(df$dailyPM, lag.max = 21)

# For the DLNM model, the first step is to 
# construct a crossbasis for each lagged variable 
# the crossbasis is a specially organized matrix 
# When we construct the crossbasis, we define 
# the function form of the relationship of the lags
# as well as the dose-response relationship 

# 4b Construct crossbasis for exposure

cb.lin.pm <- crossbasis(
  # the column with the exposure 
  df$dailyPM,         
  # The number of lags
  lag=7,    
  # the functional form of the dose-response curve
  argvar=list(fun="lin"), 
  # the functional form of the lags
  arglag=list(fun="ns", df = 6)) 


# The warning refers to this update:
# "An important change [to crossbasis() ]
# is that the knots for the spline functions or 
# cut-offs for strata in the lag dimension are now
# placed at equally-spaced percentiles if not specified, 
# differently from the default in previous versions."


# 4c Check crossbasis
# I always double check that my crossbasis makes sense 
# ie - correct number of lags, etc
# remember, summary() is a generic function, 
# and recognizes that cb.lin.pm is a crossbasis object, 
# not a model

summary(cb.lin.pm)

# 4d Construct crossbasis for temperature

cb.ns.temp <- crossbasis(
  # the column with the exposure 
  df$dailyTemp,         
  # The number of lags
  lag=7,    
  # the functional form of the dose-response curve
  argvar=list(fun="ns", df = 4), 
  # the functional form of the lags
  arglag=list(fun="ns", df = 6)) 

# 4e Construct crossbasis for relative humidity

cb.ns.rh <- crossbasis(
  # the column with the exposure 
  df$dailyRH,         
  # The number of lags
  lag=7,    
  # the functional form of the dose-response curve
  argvar=list(fun="ns", df = 4), 
  # the functional form of the lags
  arglag=list(fun="ns", df = 6)) 

# 4f Distributed Lag Model 
# note here that we are using the same 
# glm model as with the time-series analysis 
# the only thing that has changed is 
# that we are using the crossbasis terms in the formula

mod.lin <- glm(dailyCVDin ~   # outcome
                 cb.lin.pm +  # lagged, linear term for exposure
                 cb.ns.temp + # lagged, nonlinear term for Temp
                 cb.ns.rh +   # lagged, nonlinear term for RH
                 DayofWeek +  # Same-day categorical variable for DoW
                 ns(Date, df = 4*yr_num), #Same-day nonlinear term for secular trend
               family = "quasipoisson",   # distribution family
               data = df)

# 4g Extract associations 
# we can use crosspred() to estimate the association between 
# a variable and the outcome 
# at variaous values and lags
# crosspred() acts a similar role as
# pred() in previous labs

pred.lin.pm <- crosspred(
  # the exposure crossbasis
  cb.lin.pm,   
  # the model
  mod.lin, 
  # compute the estimated association for 
  # each integer value of PM2.5
  # between 0 and 55
  at = 0:55, 
  # estimates association along 
  # lags in increments of 0.2 
  # using the natural splines 
  # to interpolate between days
  bylag = 0.2,  
  # also compute cumulative association
  cumul=TRUE)

# 4h Plot results 

# 4h.i Association at each lag

plot(pred.lin.pm, 
     var=10,    # units of change of independent variable
     col="red", 
     ylab="RR", 
     ci.arg=list(density=15,lwd=2),
     main=expression("Association with a 10-unit increase in PM"[2.5]))

# Class Question ###
# How do we interpret this plot? 
####


# 4h.ii Cumulative association

plot(pred.lin.pm, "slices", var=10, col="orange", cumul=TRUE, ylab="Cumulative RR",
     main=expression("Cumulative association with a 10-unit increase in PM"[2.5]))

####***********************
#### 5: Nonlinear DLNM ####
####***********************
# We will create a new crossbasis for pm2.5
# and use the previous crossbases for temp and rh 
# since they were already non-linear

# 5a Construct crossbasis for exposure

cb.ns.pm <- crossbasis(
  # the column with the exposure 
  df$dailyPM,         
  # The number of lags
  lag = 7,    
  # the functional form of the dose-response curve
  argvar = list(fun="ns", df = 2), 
  # the functional form of the lags
  arglag = list(fun="ns", df = 4)) 

# 5b Check crossbasis

summary(cb.ns.pm)

# 5c Distributed Lag Model 

mod.ns <- glm(dailyCVDin ~   # outcome
                cb.ns.pm +    # lagged, nonlinear term for exposure
                cb.ns.temp +  # lagged, nonlinear term for Temp
                cb.ns.rh +    # lagged, nonlinear term for RH
                DayofWeek + # Same-day categorical variable for DoW
                ns(Date, df = 4*yr_num),   #Same-day nonlinear term for secular trend
              family = "quasipoisson",   # distribution family
              data = df)

# 5d Extract predictions 

pred.ns.pm <- crosspred(
  # the exposure crossbasis
  cb.ns.pm,   
  # the model
  mod.ns, 
  # compute the estimated association for 
  # each integer value of PM2.5
  # between 0 and 55
  at = 0:55, 
  # estimates association along 
  # lags in increments of 0.2 
  # using the natural splines 
  # to interpolate between days
  bylag = 0.2,  
  # center at average PM
  cen = 0,
  # also compute cumulative associations
  cumul=TRUE)

# 5e 3d Plot 

plot(pred.ns.pm, 
     xlab="\nPM2.5", zlab="\nRR", ylab="\nLag", 
     theta=40, phi=30, lphi=30,
     main=expression("3D graph of PM"[2.5]*" Association"))

# we can adjust the angle of the plot to see different perspectives 
plot(pred.ns.pm, 
     xlab="\nPM2.5", zlab="\nRR", ylab="\nLag", 
     theta=220, phi=30, lphi=30,
     main=expression("3D graph of PM"[2.5]*" Association"))

plot(pred.ns.pm, 
     xlab="\nPM2.5", zlab="\nRR", ylab="\nLag",
     theta=100, phi=30, lphi=30,
     main=expression("3D graph of PM"[2.5]*" Association"))

# 5f Slice Plot
# we can visual the association at specific lags 

# 5f.i Set colors for lags
lag.color.list <- c("firebrick", "deeppink", "mediumorchid1", "mediumpurple","purple3", "royalblue1", "darkblue")

# 5f.ii Plot dose-response curve at Lag 0
plot(pred.ns.pm, "slices",
     lag = 0,
     ci="n", 
     col=lag.color.list[1],
     ylim=c(0.95,1.10),
     lwd=1.5,
     main="Dose-Response Curves for Different Lags", 
     xlab = expression("PM"[2.5]), 
     ylab = "Rate Ratio of CVD Hospitalizations")

# 5f.iii Plot dose-response curve for Lags 1-6

for(i in 1:6) lines(pred.ns.pm, "slices", lag=c(1:6)[i], col=lag.color.list[i+1], lwd=1.5)

# 5.f.iv Add legend

legend("topleft",paste("Lag =",c(0:6)), col=lag.color.list, lwd=1.5)

# 5.f.v Add rug

rug(df$dailyPM)

####***********************************************************
#### FootNote: Examining Different Lag DF for Linear Model ####
####***********************************************************

# A Create confounder crossbases
# A1 Construct crossbasis for temperature

cb.ns.temp <- crossbasis(df$dailyTemp, lag=7, argvar=list(fun="ns", df = 4),  arglag=list(fun="ns", df = 6)) 

# A2 Construct crossbasis for relative humidity

cb.ns.rh <- crossbasis(df$dailyRH, lag=7, argvar=list(fun="ns", df = 4), arglag=list(fun="ns", df = 6)) 

# B: Construct model with 6 df for lag 
# B1 Construct crossbasis for exposure

cb.6df.pm <- crossbasis(df$dailyPM, lag=7, argvar=list(fun="lin"), arglag=list(fun="ns", df = 6)) 

# B2 Distributed Lag Model 

mod.6df <- glm(dailyCVDin ~ cb.6df.pm + cb.ns.temp + cb.ns.rh + DayofWeek + ns(Date, df = 4*yr_num),  
               family = "quasipoisson", data = df)

# B3 Extract associations 

pred.6df.pm <- crosspred(cb.6df.pm, mod.6df, at = 0:55,  bylag = 0.2,  cumul=TRUE)

# C: Construct model with 5 df for lag 
# C1 Construct crossbasis for exposure

cb.5df.pm <- crossbasis(df$dailyPM, lag=7, argvar=list(fun="lin"), arglag=list(fun="ns", df = 5)) 

# C2 Distributed Lag Model 

mod.5df <- glm(dailyCVDin ~ cb.5df.pm + cb.ns.temp + cb.ns.rh + DayofWeek + ns(Date, df = 4*yr_num),  
               family = "quasipoisson", data = df)

# C3 Extract associations 

pred.5df.pm <- crosspred(cb.5df.pm, mod.5df, at = 0:55,  bylag = 0.2,  cumul=TRUE)

# C: Construct model with 5 df for lag 
# C1 Construct crossbasis for exposure

cb.4df.pm <- crossbasis(df$dailyPM, lag=7, argvar=list(fun="lin"), arglag=list(fun="ns", df = 4)) 

# C2 Distributed Lag Model 

mod.4df <- glm(dailyCVDin ~ cb.4df.pm + cb.ns.temp + cb.ns.rh + DayofWeek + ns(Date, df = 4*yr_num),  
               family = "quasipoisson", data = df)

# C3 Extract associations 

pred.4df.pm <- crosspred(cb.4df.pm, mod.4df, at = 0:55,  bylag = 0.2,  cumul=TRUE)


# E: Construct model with 3 df for lag 
# E1 Construct crossbasis for exposure

cb.3df.pm <- crossbasis(df$dailyPM, lag=7, argvar=list(fun="lin"), arglag=list(fun="ns", df = 3)) 

# E2 Distributed Lag Model 

mod.3df <- glm(dailyCVDin ~ cb.3df.pm + cb.ns.temp + cb.ns.rh + DayofWeek + ns(Date, df = 4*yr_num),  
               family = "quasipoisson", data = df)

# E3 Extract associations 

pred.3df.pm <- crosspred(cb.3df.pm, mod.3df, at = 0:55,  bylag = 0.2,  cumul=TRUE)

# F Plot the 4 models 
#pdf(paste0(DataPath, "Df_of_Lag_for_linear_PM.pdf"))

par(mfrow = c(2,2))

plot(pred.6df.pm, var=10, col="red", ylab="RR", 
     ci.arg=list(density=15,lwd=2),
     main="6df for Lag")

plot(pred.5df.pm, var=10, col="red", ylab="RR", 
     ci.arg=list(density=15,lwd=2),
     main="5df for Lag")

plot(pred.4df.pm, var=10, col="red", ylab="RR", 
     ci.arg=list(density=15,lwd=2),
     main="4df for Lag")

plot(pred.3df.pm, var=10, col="red", ylab="RR", 
     ci.arg=list(density=15,lwd=2),
     main="3df for Lag")

#dev.off()
