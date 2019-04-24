# Time Series Lab Session
# Marianthi-Anna Kioumourtzoglou
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# Feburary 16, 2019
# Session 5: CaseCrossover 
##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Load NY Daily Data
# 2: Time Series of Air Pollution
# 3: Time Series of Cardiovascular Hospitalization 
# 4: CVD by Day of the Week 
# 5: Poisson Regression - Lag 1 
# 6: Poisson Regression - 3 Day Ave
# Footnote: Poisson Regression - Nonlinear PM 
# Footnote: Choosing df for Temp and RH 
# Footnote: Manually Plotting Time Series with Smoothing Splines 
# Footnote: Exploring Different Knots for Smoothing Splines
########################
#### 0: Preparation ####
# 0a Install packages 

# the lubridate package provides some nice functions for dealing with dates 
#install.packages("lubridate")
# 0b Load packages

library(readr)
library(dplyr) 
library(lubridate)
library(splines)
library(mgcv)
library(ggplot2)



################################
#### 1: Load NY Daily Data #####
#### 
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

head(df)
tail(df)
hist(df$dailyPM)
hist(df$dailyCVDin)
dim(df)
summary(df)

# 1c Convert date column to datetime format 
# lubridate has a nice function, parse_date_time() 
# to convert numeric or charater values to datetime format 
# all you have to do is specify the order of the 
# date information 
# ex: "dmy HM" is 
# day, month, year, hour, minute. 
 
# 1d Extract the day of the week
# lubridate has easy functions to extract the day of the week 
# the arguement label=TRUE tells wday to provide 
# the day of the week as a word, not a number 

df <- df %>% mutate(Date = parse_date_time(dateD, "mdy"))

# rearrange the rows according to descending date

df <- df %>% arrange(Date)

# day of the week

df <- df %>% mutate(DayofWeek = as.character(wday(Date, label = TRUE)))
head(df$DayofWeek)

#########################################
#### 2: Time Series of Air Pollution ####
####
# Two common ways to plot time series is 
# line plot 
# smoothed plot 

# 2a Line plot

ggplot(df, aes(Date, dailyPM)) + 
  geom_line(color = "darkorange3")


####### Class Question #### 
# Do we see any trends? 
###########################

# 2b Smoothing plot
# We can create a smoothed line across the dates 
# using penalized splines
# with the geom_smooth() function

# 2b.i Compute number of years:
yr_num <- length(unique(year(df$Date)))

ggplot(df, aes(Date, dailyPM)) + 
  geom_line(color = "darkorange2") +
  geom_smooth(color = "darkorange3", method = lm, formula = y ~ ns(x, 4 * yr_num)) +
  geom_smooth(color = "darkorange4")
  
  
##########################################################
#### 3: Time Series of Cardiovascular Hospitalization ####
####
# 3a Line plot

ggplot(df, aes(Date, dailyCVDin)) + 
  geom_line(color = "mediumpurple3")

####### Class Question #### 
# Do we see any trends? 
###########################

# 3b Smoothing plot

ggplot(df, aes(Date, dailyCVDin)) + 
  geom_smooth(color = "mediumpurple3", method = lm, formula = y ~ ns(x, 4 * yr_num)) +
  geom_smooth(color = "mediumpurple4")

###################################
#### 4: CVD by Day of the Week ####
####
# 4a Compute average number of hospitalizations per day of week

# 4a.i Order the days of the week

df$DayofWeek <- factor(df$DayofWeek, level = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"))

# 4a.ii Compute average number of hospitalizations

df.dow.cvd <- df %>% group_by(DayofWeek) %>% 
              summarize(meanCVD = mean(dailyCVDin))

df.dow.pm <- df %>% group_by(DayofWeek) %>% 
              summarize(meanPM = mean(dailyPM))

# 4b Plot 

ggplot(df.dow.cvd, aes(DayofWeek, meanCVD)) +  
  geom_bar(stat = "identity")

ggplot(df.dow.pm, aes(DayofWeek, meanPM)) +  
  geom_bar(stat = "identity")

########################################
#### 5: Poisson Regression - Lag 1  ####
####
# for this model we will not consider lagged exposure. 
# we will initially just consider same-day exposure 

# for this model we will assume linearity in the PM-CVD relationship 
# and use natural spline terms for temperature and RH 
# with 4 and 3 df. 
# For secular and seasonal trends, 
# we will use a natural spline with 4 df per year. 
# see footnotes for details. 

# 5a.ii Model 

mod.ts.p <- glm(dailyCVDin ~              # outcome
                dailyPM +                 # exposure
                ns(dailyTemp, df = 4) +   # nonlinear term for Temp
                ns(dailyRH, df = 3) +     # nonlinear term for RH
                DayofWeek +               # categorical Variable for DoW
                ns(Date, df = 4*yr_num),  # nonlinear term for secular trend
                family = "poisson"      , # distribution family
                data = df)


# 5b Summary 

summary(mod.ts.p)

# 5c Does data follow a Poisson distribution?
 
# In the poisson distribution, the variance of the data 
# is equal to the mean of the data 
# For our data, we may have overdispersion 
# where the variance is greater than the mean 
# or underdispersion, where it is lower

mean(df$dailyCVDin)
var(df$dailyCVDin)

# 5d Model -- QuasiPoisson

mod.ts.qp <- glm(dailyCVDin ~ # outcome
                       dailyPM +                # exposure
                       ns(dailyTemp, df = 4) +  # nonlinear term for Temp
                       ns(dailyRH, df = 3) +    # nonlinear term for RH
                       DayofWeek +              # categorical Variable for DoW
                       ns(Date, df = 4*yr_num), # nonlinear term for secular trend
                       family = "quasipoisson", # distribution family
                       data = df)


# 5e Summary 

summary(mod.ts.qp)

# compare model coefficients 

summary(mod.ts.p)$coef[2,]
summary(mod.ts.qp)$coef[2,]

# extract dispersion parameter 

summary(mod.ts.qp)$dispersion

# 5f Compute Rate Ratio for 1-unit increase in PM2.5

# 5f.i Create coefficient matrix

coeff.mat <- summary(mod.ts.qp)$coef

# 5f.ii Extract the coefficient and its standard error 

coeff_pm <- coeff.mat[2,1]
se_pm    <- coeff.mat[2,2]

# 5f.iii Exponentiate coefficient & CI
# Since the link function for the poisson model is log()
# we exponentiate the coefficient to compute rate ratio 

beta.pm <- exp(coeff_pm)
lci     <- exp(coeff_pm - 1.96*se_pm)
uci     <- exp(coeff_pm + 1.96*se_pm)

# 5f.iv Present results 

RR_pm <- paste(round(beta.pm, 4)," (95%CI: ", round(lci,4), ", ", round(uci,4),")", sep = "")
RR_pm

# 5g Compute Percent Change for 1-unit increase in PM2.5

perc.change <- 100*(exp(coeff_pm) - 1)
perc.lci    <- 100*(exp((coeff_pm - 1.96*se_pm)) - 1)
perc.uci    <- 100*(exp((coeff_pm + 1.96*se_pm)) - 1)

PC_pm.1 <- paste(round(perc.change, 2),"% (95%CI: ", round(perc.lci,2), "%, ", round(perc.uci,2),"%)", sep = "")
PC_pm.1

# 5h Compute Percent Change for 10-unit increase in PM2.5

perc.change <- 100*(exp(10*coeff_pm) - 1)
perc.lci    <- 100*(exp(10*(coeff_pm - 1.96*se_pm)) - 1)
perc.uci    <- 100*(exp(10*(coeff_pm + 1.96*se_pm)) - 1)

PC_pm.10 <- paste(round(perc.change, 2),"% (95%CI: ", round(perc.lci,2), "%, ", round(perc.uci,2),"%)", sep = "")
PC_pm.10

###########################################
#### 6: Poisson Regression - 3 Day Ave ####
####

# 6a Compute 3-day moving averages 
# arrange() orders the rows of our data based on the Date column 
# lag() refers to the value of the row above the current row
# lead() refers to the value of the row below the current row

# put rows in chronological order
df <- df %>% arrange(Date) %>%   
  # compute 3-day averages for pm, temp, and rh
  # for each day, compute the mean of the exposure of that day and the previous two days 
  mutate(pm_ave3d = (dailyPM + lag(dailyPM) + lag(dailyPM,2))/3,
         temp_ave3d = (dailyTemp + lag(dailyTemp) + lag(dailyTemp,2))/3,
         rh_ave3d = (dailyRH + lag(dailyRH) + lag(dailyRH,2))/3)

# 6b QuasiPoisson model with 3-day averages

mod.ts.3day.qp <- glm(dailyCVDin ~          # outcome
                   pm_ave3d +               # exposure
                   ns(temp_ave3d, df = 4) + # nonlinear term for Temp
                   ns(rh_ave3d, df = 3) +   # nonlinear term for RH
                   DayofWeek +              # categorical Variable for DoW
                   ns(Date, df = 4*yr_num), # nonlinear term for secular trend
                   family = "quasipoisson", # distribution family
                   data = df)

# 6c Summary 

# 6c.i overall summary 
summary(mod.ts.3day.qp)

# 6c.ii extract dispersion term 
summary(mod.ts.3day.qp)$dispersion

# 6c.iii extract coefficients for pm 

summary(mod.ts.3day.qp)$coef[2,]

# 6c.iv extract RR and 95% CI 

# get coef and se 

coef.pm <- summary(mod.ts.3day.qp)$coef[2,1]
se.pm <- summary(mod.ts.3day.qp)$coef[2,2]

# exponentiate coef and 95% CI 

beta.pm <- exp(coeff_pm)
lci     <- exp(coeff_pm - 1.96*se_pm)
uci     <- exp(coeff_pm + 1.96*se_pm)

# present RR and 95% CI 

RR_pm <- paste0(round(beta.pm, 4) ,
                " (95% CI: ", 
                round(lci, 4), 
                ", ",
                round(uci, 4),
                ")")
RR_pm

#################################
#### Footnote: NonLinear PM  ####
####

# A QuasiPoisson model with 3-day averages

mod.ts.3day.qp.psp <- gam(dailyCVDin ~               # outcome
                        s(pm_ave3d) +            # exposure
                        ns(temp_ave3d, df = 4) + # nonlinear term for Temp
                        ns(rh_ave3d, df = 3) +   # nonlinear term for RH
                        DayofWeek +              # categorical Variable for DoW
                        ns(Date, df = 4*yr_num), # nonlinear term for secular trend
                      family = "quasipoisson",   # distribution family
                      data = df)

# B Plot expected hospitalizations 

# B.i Create predictions for 4 degrees of freedom

predtemp <- predict(mod.ts.3day.qp.psp, se.fit = TRUE, type = "terms" )

# B.ii convert to dataframe 

predtemp <- as.data.frame(predtemp)

# B.iii Combine predictions and standard errors

predtemp <- predtemp %>% 
  rename( pred.fit = fit.s.pm_ave3d.,
          se = se.fit.s.pm_ave3d.)

# B.iv Compute 95% confidence intervals 

predtemp <- predtemp %>% 
  mutate( pred = exp(pred.fit),
          lci = exp(pred.fit - 1.96*se),
          uci = exp(pred.fit + 1.96*se))

# B.v Keep only variables we need

predtemp <- predtemp %>% select(pred, se, lci, uci) 

# B.vi Combine with data 

df.pred <- df[3:nrow(df),] %>% bind_cols(predtemp)

# B.vii Uncenter the data - multiply by mean dailyCVDin

# the predicted values should be multiplied
# because the poisson model is multiplicative at the count scale 
# ie - a prediction of 0.7 means that, 
# this model predicts that that observation
# given its observed PM, and median value of other variables
# has a log(hospitalization rate) that is 0.7 greater 
# than the average log(hospitalization rate)
# log(PredHRate) = 0.7 + log(AveHRate)
# PredHRate = exp(0.7 + log(AveHRate))
# PredHRate = exp(0.7) * exp(log(AveHRate))
# PredHRate = exp(0.7) * aveHRate

df.pred <- df.pred %>% mutate(predPM = pred * mean(dailyCVDin),
                              lciPM = lci * mean(dailyCVDin),
                              uciPM = uci * mean(dailyCVDin))

# B.viii Plot 

ggplot(df.pred, aes(pm_ave3d)) + 
  geom_line(aes(y = predPM), color = "darkorange3") + 
  geom_line(aes(y = lciPM), color = "grey", alpha = 0.5) + 
  geom_line(aes(y = uciPM), color = "grey", alpha = 0.5) + 
  xlab(expression("3-Day Moving Average PM"[2.5])) +
  ylab("Expected Number of Hospitalizations")  


###############################################
#### Footnote: Choosing df for Temp and RH ####
####
# we will choose our df based on the relationship between 
# the confounders and the exposure of interest 
# we will chose the term with the lowest AIC 

# A Create models for temperature
# Here we will use a natural spline term for dailyRH 

mod.lin  <- lm(dailyPM ~ dailyTemp + ns(dailyRH, df = 2) + DayofWeek, data = df)
mod.ns.2 <- lm(dailyPM ~ ns(dailyTemp, df = 2) + ns(dailyRH, df = 2)+ DayofWeek, data = df)
mod.ns.3 <- lm(dailyPM ~ ns(dailyTemp, df = 3) + ns(dailyRH, df = 2) + DayofWeek, data = df)
mod.ns.4 <- lm(dailyPM ~ ns(dailyTemp, df = 4) + ns(dailyRH, df = 2) + DayofWeek, data = df)

# B Create AIC table for temperature

models_aic        <- data.frame( c("linear", "2 df", "3 df", "4 df"),
                                 c(AIC(mod.lin), AIC(mod.ns.2), AIC(mod.ns.3), AIC(mod.ns.4)), 
                                 stringsAsFactors = FALSE)
names(models_aic) <- c("ModelName", "AIC")
models_aic$ModelName[which(models_aic$AIC==min(models_aic$AIC))]

# 4 df is selected for the natural spline term for Temperature

# C Create models for RH
# Here we will use a natural spline term for dailyRH 

mod.lin  <- lm(dailyPM ~ ns(dailyTemp, df = 4) + dailyRH + DayofWeek, data = df)
mod.ns.2 <- lm(dailyPM ~ ns(dailyTemp, df = 4) + ns(dailyRH, df = 2) + DayofWeek, data = df)
mod.ns.3 <- lm(dailyPM ~ ns(dailyTemp, df = 4) + ns(dailyRH, df = 3) + DayofWeek, data = df)
mod.ns.4 <- lm(dailyPM ~ ns(dailyTemp, df = 4) + ns(dailyRH, df = 4) + DayofWeek, data = df)

# D Create AIC table for RH
models_aic        <- data.frame( c("linear", "2 df", "3 df", "4 df"),
                                 c(AIC(mod.lin), AIC(mod.ns.2), AIC(mod.ns.3), AIC(mod.ns.4)), 
                                 stringsAsFactors = FALSE)
names(models_aic) <- c("ModelName", "AIC")
models_aic$ModelName[which(models_aic$AIC==min(models_aic$AIC))]

# 3 df is selected for the natural spline term for Relative Humidity


########################################################################
#### Footnote: Manually Plotting Time Series with Smoothing Splines ####

# A Construct smoothing model 
# 6 degrees of freedom per year

mod.temp.smooth <- lm(dailyTemp ~ ns( Date, df = 48), data= df)

# B Create predictions for 4 degrees of freedom

predtemp <- predict(mod.temp.smooth, se.fit = TRUE, type = "terms" )

# C convert to dataframe 

predtemp <- as.data.frame(predtemp)
names(predtemp)

# D Combine predictions and standard errors

predtemp <- predtemp %>% 
  rename( pred = ns.Date..df...48.,
          se = ns.Date..df...48..1)

# E Compute 95% confidence intervals 

predtemp <- predtemp %>% 
  mutate( lci = pred - 1.96*se,
          uci = pred + 1.96*se)

# F Keep only variables we need

predtemp <- predtemp %>% select(pred, se, lci, uci) %>% 
  mutate(Model = "48 df")

# G Combine with data 

df.pred <- df %>% bind_cols(predtemp)

# F Uncenter the data - add mean aveBMI

df.pred <- df.pred %>% mutate(predTemp = pred + mean(dailyTemp),
                              lciTemp = lci + mean(dailyTemp),
                              uciTemp = uci + mean(dailyTemp))

# G Plot 

ggplot(df.pred, aes(Date)) + 
  geom_line(aes(y = predTemp), color = "cadetblue3") + 
  geom_line(aes(y = lciTemp), color = "grey", alpha = 0.5) + 
  geom_line(aes(y = uciTemp), color = "grey", alpha = 0.5) + 
  xlab("Date") + 
  ylab("Predicted Temperature") 

###################################################################
#### Footnote: Exploring Different Knots for Smoothing Splines ####

ggplot(df.pred, aes(Date)) + 
  geom_point(aes(y = dailyTemp), color = "darkblue")+
  geom_line(aes(y = predTemp), color = "cadetblue3") + 
  geom_line(aes(y = lciTemp), color = "grey", alpha = 0.5) + 
  geom_line(aes(y = uciTemp), color = "grey", alpha = 0.5) + 
  geom_smooth(aes(y=dailyTemp),method = gam, formula = y ~ ns(x, df=48), color = "purple") +
  geom_smooth(aes(y=dailyTemp),method = gam, formula = y ~ s(x, k=20), color = "firebrick2") +
  xlab("Date") + 
  ylab("Predicted Temperature") 

ggplot(df.pred, aes(Date)) + 
  geom_point(aes(y = dailyTemp), color = "darkblue")+
  geom_smooth(aes(y=dailyTemp),method = gam, formula = y ~ s(x, k=10), color = "firebrick4") +
  geom_smooth(aes(y=dailyTemp),method = gam, formula = y ~ s(x, k=15), color = "firebrick3") +
  geom_smooth(aes(y=dailyTemp),method = gam, formula = y ~ s(x, k=20), color = "firebrick2") +
  #geom_smooth(aes(y=dailyTemp),method = gam, formula = y ~ s(x, k=25), color = "firebrick1") +
  geom_smooth(aes(y=dailyTemp),method = gam, formula = y ~ s(x, k=17), color = "indianred3") +
  xlab("Date") + 
  ylab("Predicted Temperature") 
