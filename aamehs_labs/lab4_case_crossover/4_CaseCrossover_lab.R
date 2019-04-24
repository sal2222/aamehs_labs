# CaseCrossover Lab Session
# Marianthi-Anna Kioumourtzoglou
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# Feburary 4, 2019
# Session 4: CaseCrossover 
##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Load County Data 
# 2: Logistic Model 
# 3: Load NYC Daily Data
# 4: CaseCrossover
# 5: Nonlinearity in CaseCrossover 
# Footnote: choosing df for temp and rh 
# Footnote: clogit vs coxph
# Footnote: Plotting Natural Splines
########################
#### 0: Preparation ####
########################
# 0a Install packages 

# the survival package provides us with the model commands for 
# conditional logistic models 
#install.packages("survival")

# 0b Load packages

library(readr)
library(dplyr) 
library(survival)
library(splines)
library(pspline)
library(ggplot2)


#############################
#### 1: Load County Data ####
#############################
# 1a Load data 

df <- read_csv(".csv")


# 1b Remove the variables we will not need 

df <- df %>% select(-AREA, -PERIMETER, -geometry)

# 1c Keep only complete cases 

df <- df %>% filter(complete.cases(df))

###########################
#### 2: Logistic Model ####
###########################
# 2a Create a dichotomous variable 
# for whether average county bmi is greate than average

df <- df %>% mutate(aveObese = if_else(aveBMI > 29.6, 1, 0))

# 2b Create a logistic model 

mod.log <- glm(aveObese ~ avePM.idw + maleunemp  + femaleunemp + ltHS + 
               medhinc + medhval + per.black + per.latinx + per.asnam + climate_region,
               family = binomial(link=logit),
               data = df)

# 2c Model Summary 

summary(mod.log)

# 2d Compute odds ratio and 95% CI 
# for a 1-unit increase in pm2.5 

# 2d.i Create coefficient matrix

coeff.mat <- summary(mod.log)$coef

# 2d.ii Extract the coefficient and standard error 

coeff_pm <- coeff.mat[2,1]
se_pm    <- coeff.mat[2,3]

# 2d.iii Exponentiate coefficient & CI
# since the logistic model estimates the log of the odds
# we can compute the odds ratio by exponentiating the coefficients 
# in R exp(a) is the same as 'e to the power of a'

beta.pm <- exp(coeff_pm)
lci     <- exp(coeff_pm - 1.96*se_pm)
uci     <- exp(coeff_pm + 1.96*se_pm)

# 2d.iv Present results 

OR_pm <- paste(round(beta.pm, 2)," (95%CI: ", round(lci,2), ", ", round(uci,2),")", sep="")
OR_pm

# 2e Compute percent change in risk and 95% CI 
# when an event or disease is rare, 
# the odds ratio approximates the risk ratio 

perc.change <- 100*(exp(coeff_pm)-1)
perc.lci    <- 100*(exp((coeff_pm - 1.96*se_pm))-1)
perc.uci    <- 100*(exp((coeff_pm + 1.96*se_pm))-1)

PC_pm <- paste(round(perc.change, 2),"% (95%CI: ", round(perc.lci,2), "%, ", round(perc.uci,2),"%)", sep="")
PC_pm

# 2.f Clean the environment

rm(df, coeff.mat, beta.pm, coeff_pm, lci, mod.log, OR_pm,
   perc.change, perc.lci, perc.uci, PC_pm, se_pm, uci)

################################
#### 3: Load NYC Daily Data ####
################################
# 3a Readin data 

df <- read_csv("nyc_daily_processed.csv")

# 3b Keep only complete cases 

df <- df %>% filter(complete.cases(df))

# 3c Look at data structure 
# Data Dictionary 
# all rows are either case days or control days 
# DayID - date of the corresponding Case Day 
# DayName - whether it is a case, 
# or, if a control, how many weeks before or after case
# Case - 1 if a case, 0 if a control 
# TimetoEvent - 1 if a case, 2 if a control 
# DailyPM - population-weighted average of AQS estimates across NYC 
# DailyTemp - population-weight average of temperature in Celcius from NLDAS 
# Daily RH - like DailyTemp, but for relative humidity
# DailyCVDin - number of in-patient hospitalizations for cardiovascular disease

df <- df %>% arrange(DayID)

head(df)
tail(df)
hist(df$dailyPM)

dim(df)
summary(df)
table(df$Case)
length(unique(df$DayDate))

###########################
# Simple time-series for PM trend
###########################
Cases <- df %>% filter(Case ==1)

ggplot(df, aes(DayDate, dailyPM))+ 
  geom_smooth(color = "cyan")


###########################
#### 4: CaseCrossover  ####
###########################
# we use the clogit command to create a conditional logistic model 
# we will start with the assumption that dailyPM is linear
# we already expect that temperature and relative humidity will 
# have a non-linear relationship to the outcome 
# so we include a natural spline term for temperature and RH 
# The degrees of freedom were chosen based on AIC
# see footnote for more details

# 4a Model

mod.clogit.lin <- clogit(Case ~ dailyPM + ns(dailyTemp, df=3) + ns(dailyRH, df=4) + 
               strata(DayID), # each case day is a strata
               weights = dailyCVDin, #number of events in each day
               method = "efron", # the method tells the model how to deal with ties
               df) 

# 4b Model Summary 

summary(mod.clogit.lin)

# 4c Extract odds ratio 

# 4c.i Create coefficient matrix

coeff.mat <- summary(mod.clogit.lin)$coef

# 4c.ii Extract the coefficient and standard error 

coeff_pm <- coeff.mat[1,1]
se_pm <- coeff.mat[1,3]

# 4c.iii Exponentiate coefficient & CI
# since the conditional logistic model involves taking the log of the odds
# we can compute the odds ratio by exponentiating the coefficients 
# in R exp(a) is the same as 'e to the power of a'

beta.pm <- exp(coeff_pm)
lci     <- exp(coeff_pm - 1.96*se_pm)
uci     <- exp(coeff_pm + 1.96*se_pm)

# 4c.iv Present results 

OR_pm <- paste(round(beta.pm, 4)," (95%CI: ", round(lci,4), ", ", round(uci,4),")", sep="")
OR_pm

# 4d Compute Odds Ratio for 10-unit increase
# we might also be interested in how the odds vary 
# for larger increments of pm2.5 

# 4d.i compute beta and CI 

beta.pm.10 <- exp(10 * coeff_pm)
lci.10     <- exp(10 * (coeff_pm - 1.96*se_pm))
uci.10     <- exp(10 * (coeff_pm + 1.96*se_pm))

# 4c.iv Present results 

OR_pm.10 <- paste(round(beta.pm.10, 4)," (95%CI: ", round(lci.10,4), ", ", round(uci.10,4),")", sep="")
OR_pm.10


# 4d Plot model 
# 4d.i Predict odds and standard error

ptemp <- predict(mod.clogit.lin, type = "terms", se = TRUE)
ptemp <- as.data.frame(ptemp) 

colnames(ptemp) <- c("term.PM","termTemp.ns.3", "termRH.ns.3", "se.PM","seTemp.ns.3", "seRH.ns.3")
df.pred         <- df %>% bind_cols(ptemp)

# 4d.ii Compute confidence intervals 

df.pred <- df.pred %>%
           mutate(lci = term.PM - 1.96 * se.PM , 
                  uci = term.PM + 1.96 * se.PM )

# 4d.iii Exponentiate 
df.pred <- df.pred %>%
           mutate(termPM = exp(term.PM),
           lci = exp(lci), 
           uci = exp(uci) )

# 4d.iv Plot

ggplot(df.pred, aes(dailyPM)) + 
  geom_line(aes(y = termPM), color = "lightsalmon2", size = 1) +
  geom_line(aes(y = lci), color = "grey", linetype = "dashed") + 
  geom_line(aes(y = uci), color = "grey", linetype = "dashed") +  
  geom_rug(col = "gray60") +
  xlab(expression("Daily PM"[2.5]))+
  ylab("Odds of CVD admission")


################################################
#### Footnote: Choosing df for Confounders  ####
################################################
# we will first identify the best-fitting degrees of freedom for 
# the temperature-cvd relationship 
# and then for rh 
# Process for df selection 
# 1: eliminate any df that lead to models that are too wiggly to be plausible 
# 2: choose the df that leads to the lowest aic. 
# 1: Temperature
# 1A Create the models

mod.lin <- clogit(Case ~  dailyTemp + strata(DayID), weights = dailyCVDin,
                 method = "efron", df) 
mod.ns.2 <- clogit(Case ~  ns(dailyTemp, df=2) + strata(DayID), weights = dailyCVDin,
               method = "efron", df) 
mod.ns.3 <- clogit(Case ~  ns(dailyTemp, df=3) + strata(DayID), weights = dailyCVDin,
                   method = "efron", df) 
mod.ns.4 <- clogit(Case ~  ns(dailyTemp, df=4) + strata(DayID), weights = dailyCVDin,
                   method = "efron", df) 

# 1B create predictions for each model

df.pred <- df

ptemp <- predict(mod.lin, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.lin")
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.2, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.ns.2")
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.3, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.ns.3")
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.4, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.ns.4")
df.pred <- df.pred %>% bind_cols(ptemp)

# 1C Plot 
ggplot(df.pred, aes(dailyTemp)) + 
  geom_line(aes(y = exp(termT.lin)), color = "springgreen3", size = 1) +
  geom_line(aes(y = exp(termT.ns.2)), color = "Lightsalmon", size = 1) +
  geom_line(aes(y = exp(termT.ns.3)), color = "lightsalmon2", size = 1) +
  geom_line(aes(y = exp(termT.ns.4)), color = "lightsalmon4", size = 1) +
  xlab("dailyTemp")+
  ylab("Odds of CVD admission")

# none of the terms are too wiggly, so we will base our decision just on aic 
# 1D AIC

models_aic        <- data.frame( c("linear", "2 df", "3 df", "4 df"),
                                 c(AIC(mod.lin), AIC(mod.ns.2), AIC(mod.ns.3), AIC(mod.ns.4)), 
                                 stringsAsFactors = FALSE)
names(models_aic) <- c("ModelName", "AIC")
models_aic$ModelName[which(models_aic$AIC==min(models_aic$AIC))]

# 3 df is selected for the natural spline term for meanT
# 2 repeat the process for dailyRH. 

# 2A Create the models

mod.lin <- clogit(Case ~  dailyRH + strata(DayID), weights = dailyCVDin,
                  method = "efron", df) 
mod.ns.2 <- clogit(Case ~  ns(dailyRH, df=2)+strata(DayID), weights = dailyCVDin,
                   method = "efron", df) 
mod.ns.3 <- clogit(Case ~  ns(dailyRH, df=3)+strata(DayID), weights = dailyCVDin,
                   method = "efron", df) 
mod.ns.4 <- clogit(Case ~  ns(dailyRH, df=4)+strata(DayID), weights = dailyCVDin,
                   method = "efron", df) 

# 2B create predictions for each model
df.pred <- df

ptemp <- predict(mod.lin, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.lin")
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.2, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.ns.2")
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.3, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.ns.3")
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.4, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c("termT.ns.4")
df.pred <- df.pred %>% bind_cols(ptemp)

# 2C Plot 

ggplot(df.pred, aes(dailyRH)) + 
  geom_line(aes(y = exp(termT.lin)), color = "springgreen3", size = 1) +
  geom_line(aes(y = exp(termT.ns.2)), color = "Lightsalmon", size = 1) +
  geom_line(aes(y = exp(termT.ns.3)), color = "lightsalmon2", size = 1) +
  geom_line(aes(y = exp(termT.ns.4)), color = "lightsalmon4", size = 1) +
  xlab("dailyRH")+
  ylab("Odds of CVD admission")

# the 4 df model is odd, and is likely overfitting 
# 2D AIC

models_aic        <- data.frame( c("linear", "2 df", "3 df", "4 df"),
                                 c(AIC(mod.lin), AIC(mod.ns.2), AIC(mod.ns.3), AIC(mod.ns.4)), 
                                 stringsAsFactors = FALSE)
names(models_aic) <- c("ModelName", "AIC")
models_aic$ModelName[which(models_aic$AIC==min(models_aic$AIC))]

# 4 df is selected for the natural spline term for RH

###################################
#### Footnote: coxph vs clogit ####
###################################
# create the two models
mod.coxph <- coxph(Surv(TimetoEvent, Case) ~ dailyPM + 
                     ns(dailyTemp, 3) + ns(dailyRH, 3)+
                     strata(DayID), 
                   weights = dailyCVDin,
                   df)

mod.clogit <- clogit(Case ~ dailyPM + ns(dailyTemp, 3) + ns(dailyRH, 3) + strata(DayID), 
                     weights = dailyCVDin,
                     method = "efron",
                     df)
# compare their coefficients 
summary(mod.coxph)
summary(mod.clogit)

#############################################
#### Footnote: Plotting Natural Splines  ####
#############################################
# 1A Create the models

mod.ns.2 <- clogit(Case ~  ns(dailyPM, df = 2) + ns(dailyTemp, df=3) + ns(dailyRH, df=4) + 
                     strata(DayID), weights = dailyCVDin, method = "efron", df)  
mod.ns.3 <- clogit(Case ~  ns(dailyPM, df = 3) + ns(dailyTemp, df=3) + ns(dailyRH, df=4) + 
                     strata(DayID), weights = dailyCVDin, method = "efron", df) 
mod.ns.4 <- clogit(Case ~  ns(dailyPM, df = 4) + ns(dailyTemp, df=3) + ns(dailyRH, df=4) + 
                     strata(DayID), weights = dailyCVDin, method = "efron", df)  

# 1B create predictions for each model
df.pred <- df

ptemp <- predict(mod.ns.2, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c(c("termPM.ns.2"), colnames(ptemp)[2:ncol(ptemp)])
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.3, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c(c("termPM.ns.3"), colnames(ptemp)[2:ncol(ptemp)])
df.pred <- df.pred %>% bind_cols(ptemp)

ptemp <- predict(mod.ns.4, type = "terms")
ptemp <- as.data.frame(ptemp) 
colnames(ptemp) <- c(c("termPM.ns.4"), colnames(ptemp)[2:ncol(ptemp)])
df.pred <- df.pred %>% bind_cols(ptemp)

# 1C Plot 
pdf(paste0(output.folder, "NaturalSplines_for_PM.pdf"))
ggplot(df.pred, aes(dailyPM)) + 
  geom_line(aes(y = exp(termPM.ns.2)), color = "Lightsalmon", size = 1) +
  geom_line(aes(y = exp(termPM.ns.3)), color = "lightsalmon2", size = 1) +
  geom_line(aes(y = exp(termPM.ns.4)), color = "lightsalmon4", size = 1) +
  xlab("dailyPm")+
  ylab("Odds of CVD admission")
dev.off()