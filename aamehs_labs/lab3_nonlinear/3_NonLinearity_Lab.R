# Non-Linearity Lab Session
# Marianthi-Anna Kioumourtzoglou
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# Feburary 4, 2019
# Session 3: Non-Linearity 
##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Load Data 
# 2: Quadratic Term
# 3: Piecewise Linear Spline Term 
# 4: Natural Spline Term 
# 5: Penalized Spline Term 
# Footnote: Plotting All the Models Together
# Footnote: Cubic Term 
# Footnote: Multiple Plots
########################
#### 0: Preparation ####
########################
# 0a Install packages 

# mgcv allows for quick modeling of generalized additive models 
# install.packages("mgcv")

# splines gives additional spline options
# install.packages("splines")

# 0b Load packages
library(tidyverse)
library(mgcv)
library(splines)


###############################################################
#######################
##### 1: Load Data ####
#######################
# 1a Load data 


df <- read_csv("C:/Users/slewa/Desktop/AAMEHS/aamehs_labs/lab3_nonlinear/county_bmi_pm_confounders_10.csv")

# 1b Remove the variables we will not need 

df <- df %>% select(-AREA, -PERIMETER, -geometry)

# 1c Keep only complete cases 

df <- df %>% filter(complete.cases(df))

# alternative command 
# df <- na.omit(df)

############################
##### 2: Quadratic Term ####
############################

# 2a Create model 
# we can create new terms within the model statement 
# using the I() command 

mod.quad <- lm(aveBMI ~ avePM.idw + I(avePM.idw^2) + maleunemp  + femaleunemp + ltHS + 
               medhinc + medhval + per.black + per.latinx + per.asnam + climate_region, 
               data = df)

# 2b Model Summary 

summary(mod.quad)
plot(mod.quad)

# 2c Construct predictions based on the model
# here predict() computes expected aveBMI for the observation 
# if that observation had the same value of that variable 
# and the mean value for all other variables 
# These predictions are also centered around 
# The mean aveBMI.
# But they are not scaled

predBMI.quad <- predict(mod.quad, se.fit = TRUE, type = "terms" )

# 2d Convert to dataframe 

predBMI.quad <- as.data.frame(predBMI.quad)

# 2e Combine predictions 

predBMI.quad <- predBMI.quad %>% 
  mutate( pred = fit.avePM.idw + fit.I.avePM.idw.2.)

# 2f Keep only variables we need 

predBMI.quad2 <- predBMI.quad %>% select(pred)

# 2g Combine with data 

predBMI.quad2 <- predBMI.quad2 %>% bind_cols(df)

# 2h Uncenter the expectations - add mean BMI 

predBMI.quad2 <- predBMI.quad2 %>% mutate(predBMI = pred + mean(aveBMI))

# 2i Plot

ggplot(predBMI.quad2, aes(x = avePM.idw)) + 
  geom_line(aes(y = predBMI)) + 
  xlab(expression("Average Annual PM"[2.5])) + 
  ylab("Predicted Average BMI") + 
  ylim(26, 29)

# 2k Assess model fit 

# Did adding the quadratic term improve the fit of our model? 
# we can compare the fit of the quadratic model 
# with the fit of the more simple linear model 

# 2k.i create linear model 

mod.lin <- lm(aveBMI ~ avePM.idw +  maleunemp  + femaleunemp + ltHS + 
              medhinc + medhval + per.black + per.latinx + per.asnam + climate_region, 
              data = df)

# 2k.ii Likelihood Ratio Test 
# anova() is another generic command, like summary or plot
# If anova() detects that its inputs are two nested models 
# then anova() will conduct a likelihood ratio test 

anova(mod.quad, mod.lin)

# Class Question ######################################
# Based on the LRT, did the quadratic term improve the model fit? 
#**********************
#yes


#############################################
##### 3: Piecewise Linear Spline Term ########
##############################################
# Another method is to create multiple individual linear terms 

# 3a Create spline term 
# the knot is at 12 
# for observations with pm <12, pm.lt12 is 0

df <- df %>% mutate( pm.gt12 = if_else(avePM.idw > 12, avePM.idw - 12, 0))

# 3b Create model with linear spline term 

mod.pls <- lm(aveBMI ~ avePM.idw + pm.gt12 + maleunemp  + femaleunemp + ltHS + medhinc + medhval + 
              per.black + per.latinx + per.asnam + climate_region, data = df)

summary(mod.pls)
# 3c Create predictions

predBMI.pls <- predict(mod.pls, se.fit = TRUE, type = "terms" )

# 3d Convert to dataframe 

predBMI.pls <- as.data.frame(predBMI.pls)

# 3e Combine predictions 

predBMI.pls <- predBMI.pls %>% 
               mutate( pred = fit.avePM.idw + fit.pm.gt12)

# 3f Keep only variables we need 

predBMI.pls2 <- predBMI.pls %>% select(pred)

# 3g Combine with data 

predBMI.pls2 <- predBMI.pls2 %>% bind_cols(df)

# 3h Uncenter the data - add mean BMI 

predBMI.pls2 <- predBMI.pls2 %>% mutate(predBMI = pred + mean(aveBMI))

# 3i Plot

ggplot(predBMI.pls2, aes(avePM.idw)) + 
  geom_line(aes(y = predBMI)) + 
  xlab(expression("Average Annual PM"[2.5])) + 
  ylab("Predicted Average BMI") + 
  ylim(26, 29)

#####################################
##### 4: Natural Spline Term ########
#####################################
# we can constuct natural splines within lm() 
# with the ns() command 
# we have to define the degrees of freedom 

# 4a Create model 

mod.ns.3 <- lm(aveBMI ~ ns(avePM.idw, df = 3) + maleunemp  + femaleunemp + ltHS + medhinc + medhval + 
             per.black + per.latinx + per.asnam + climate_region, data = df)

# 4b Model Summary 

summary(mod.ns.3)

# 4c construct predictions based on the model

predBMI.ns.3 <- predict(mod.ns.3, se.fit = TRUE, type = "terms" )

# 4d convert to dataframe 

predBMI.ns.3 <- as.data.frame(predBMI.ns.3)

# 4e Rename predictions and standard errors
# column has different names since it is one term 

predBMI.ns.3 <- predBMI.ns.3 %>% 
              mutate( pred = fit.ns.avePM.idw..df...3.,
              se = se.fit.ns.avePM.idw..df...3.)

# 4f Compute 95% confidence intervals 

predBMI.ns.3 <- predBMI.ns.3 %>% 
                mutate( lci = pred - 1.96*se,
                        uci = pred + 1.96*se)

# 4g Keep only variables we need

predBMI.ns.3.2 <- predBMI.ns.3 %>% select(pred, se, lci, uci)

# 4h Combine with data 

predBMI.ns.3.2 <- predBMI.ns.3.2 %>% bind_cols(df)

# 4i Uncenter data 

predBMI.ns.3.2 <- predBMI.ns.3.2 %>% mutate(predBMI = pred + mean(aveBMI),
                                lciBMI = lci + mean(aveBMI),
                                uciBMI = uci + mean(aveBMI))
# 4j Plot

ggplot(predBMI.ns.3.2, aes(avePM.idw)) + 
      geom_line(aes(y = predBMI)) + 
      geom_line(aes(y = lciBMI), color = "darkgrey") + 
      geom_line(aes(y = uciBMI), color = "darkgrey") + 
      xlab(expression("Average Annual PM"[2.5])) + 
      ylab("Predicted BMI (95% CI)") + 
      ylim(26, 29)

# 4k Compare models with different degrees of freedom 

# 4k.i Construct models 

mod.ns.2 <- lm(aveBMI ~ ns(avePM.idw, df = 2) + maleunemp  + femaleunemp + ltHS + medhinc + medhval + 
                per.black + per.latinx + per.asnam + climate_region, data = df)

mod.ns.4 <- lm(aveBMI ~ ns(avePM.idw, df = 4) + maleunemp  + femaleunemp + ltHS + medhinc + medhval + 
                per.black + per.latinx + per.asnam + climate_region, data = df)


# 4k.ii Extract AIC of each model 

aic.mod.ns.2 <- AIC(mod.ns.2)
aic.mod.ns.3 <- AIC(mod.ns.3)
aic.mod.ns.4 <- AIC(mod.ns.4)

# 4k.iii Put AIC's into a table 

models_aic        <- data.frame( c("2 df", "3 df", "4 df"),
                                 c(aic.mod.ns.2, aic.mod.ns.3, aic.mod.ns.4))
names(models_aic) <- c("ModelName", "AIC")
models_aic

# Class Question ######################################
# Based on the AIC, which model would we prefer? 
# Looking at the plots, which model looks more plausible? 
#**********************

# 4k.iv Select the model with the lowest AIC
# here, which() identifies rows that meet the criteria 
# and the criteria is having 
# AIC equal to the minimum AIC 
# among the models in the table 

models_aic$ModelName[which(models_aic$AIC == min(models_aic$AIC))]


# 4L Create natural spline with gam()

mod.ns.gam.3 <- gam(aveBMI ~ ns(avePM.idw, df = 3) + maleunemp  + femaleunemp + ltHS + medhinc + medhval + 
                 per.black + per.latinx + per.asnam + climate_region, data = df)
termplot(mod.ns.gam.3, se = TRUE)


########################################
##### 5: Penalized Spline Term  ########
########################################
# 5a Construct model
# we can constuct penalized splines within gam() 
# default is 10 knots 
# only increase if need >~8.5 edf 
# (estimated degrees of freedom, vs. user-defined)
# penalty is estimated -- the model selects the penalty 
# that leads to the highest gcv

mod.ps <- gam(aveBMI ~ s(avePM.idw) + maleunemp  + femaleunemp + ltHS + medhinc + medhval + 
                per.black + per.latinx + per.asnam + climate_region, data = df)

# 5b Model Summary 

summary(mod.ps)

# 5c Extract Penalty 
# this is the penalty estimated by the model 

mod.ps$sp

# 5d Model Plot 
# plot.gam offers a nice default plot
# a quick way to plot the change in predicted aveBMI with avePM.idw

plot(mod.ps)

# 5e Plot with ggplot 
# we can also recreate this plot with ggplot

# 5e construct predictions based on the model

predBMI.ps <- predict(mod.ps, se.fit = TRUE, type = "terms" )

# 5f convert to dataframe 

predBMI.ps <- as.data.frame(predBMI.ps)

# 5g Combine predictions and standard errors
# it just has different names since its one term 

predBMI.ps <- predBMI.ps %>% 
  mutate( pred = fit.s.avePM.idw.,
          se = se.fit.s.avePM.idw.)

# 5h Compute 95% confidence intervals 

predBMI.ps <- predBMI.ps %>% 
  mutate( lci = pred - 1.96*se,
          uci = pred + 1.96*se)

# 5i Keep only variables we need 

predBMI.ps2 <- predBMI.ps %>% select(pred, se, lci, uci)

# 5j Combine with data 

predBMI.ps2 <- predBMI.ps2 %>% bind_cols(df)

# 5k Uncenter data 

predBMI.ps2 <- predBMI.ps2 %>% mutate(predBMI = pred + mean(aveBMI),
                                      lciBMI = lci + mean(aveBMI),
                                      uciBMI = uci + mean(aveBMI))
# 5l Plot

ggplot(predBMI.ps2, aes(avePM.idw)) + 
  geom_line(aes(y = predBMI)) + 
  geom_line(aes(y = lciBMI), color = "darkgrey") + 
  geom_line(aes(y = uciBMI), color = "darkgrey") + 
  xlab(expression("Average Annual PM"[2.5])) + 
  ylab("Predicted BMI (95% CI)") + 
  ylim(26, 29)

#########################################################
##### Footnote: Plotting All the Models Together ########
#########################################################
# create a column of TermName for each dataframe of predictions 
# we will use this column to keep track of which model 
# the predictions came from 

# A Create linear model 
predBMI.lin  <- predict(mod.lin, se.fit = TRUE, type = "terms" )
predBMI.lin  <- as.data.frame(predBMI.lin)
predBMI.lin  <- predBMI.lin %>% mutate(pred = fit.avePM.idw)
predBMI.lin2 <- predBMI.lin %>% select(pred)
predBMI.lin2 <- predBMI.lin2 %>% bind_cols(df)
predBMI.lin2 <- predBMI.lin2 %>% mutate(predBMI = pred + mean(aveBMI))

# B Create a column with the model name 
# we will use this to keep track of the models once we combine them
predBMI.lin3  <- predBMI.lin2  %>% mutate(ModelName = "Linear")
predBMI.quad3 <- predBMI.quad2 %>% mutate(ModelName = "Quadratic Term")
predBMI.pls3  <- predBMI.pls2  %>% mutate(ModelName = "Piecewise Linear Spline")
predBMI.ns.3.3<- predBMI.ns.3.2%>% mutate(ModelName = "Natural Spline 3 df")
predBMI.ps3   <- predBMI.ps2   %>% mutate(ModelName = "Penalized Spline")

# C Combine the predictions 

predBMI.tot <- bind_rows(predBMI.lin3,
                         predBMI.quad3, 
                         predBMI.pls3,
                         predBMI.ns.3.3,
                         predBMI.ps3)
# D Plot!

all.models.plot <-  ggplot(predBMI.tot, aes(avePM.idw)) + 
  geom_line(aes(y = predBMI, color = ModelName)) + 
  xlab(expression("Average Annual PM"[2.5])) + 
  ylab("Predicted BMI") + 
  ylim(26, 29)

all.models.plot 

# save plot as a pdf
pdf(paste0(OutputPath, "all_models_plot.pdf"), width = 10)
all.models.plot
dev.off()

###################################
##### Footnote: Cubic Term ########
###################################
# A Create model 
# we can create new terms within the model statement 
# using the I() command 

mod.cub <- lm(aveBMI ~ avePM.idw + I(avePM.idw^2) + I(avePM.idw^3) + 
              maleunemp  + femaleunemp + ltHS + medhinc + medhval +
              per.black + per.latinx + per.asnam + climate_region, 
              data = df)

# B Model Summary 

summary(mod.cub)

# C Construct predictions based on the model

predBMI <- predict(mod.cub, se.fit = TRUE, type = "terms" )

# D convert to dataframe 

predBMI <- as.data.frame(predBMI)

# E Combine predictions and standard errors

predBMI <- predBMI %>% 
  mutate( pred = fit.avePM.idw + fit.I.avePM.idw.2.+ fit.I.avePM.idw.3.,
          se = se.fit.avePM.idw + se.fit.I.avePM.idw.2.+ se.fit.I.avePM.idw.3.)

# F Compute 95% confidence intervals 

predBMI <- predBMI %>% 
  mutate( lci = pred - 1.96*se,
          uci = pred + 1.96*se)

# G Keep only variables we need 

predBMI2 <- predBMI %>% select(pred, se, lci, uci)

# H Combine with data 

predBMI2 <- predBMI2 %>% bind_cols(df)

# I Plot

ggplot(predBMI2, aes(avePM.idw)) + 
  geom_line(aes(y = pred)) + 
  geom_line(aes(y = lci), color = "darkgrey") + 
  geom_line(aes(y = uci), color = "darkgrey") + 
  xlab(expression("Average Annual PM"[2.5])) + 
  ylab("Predicted aveBMI")

#######################################
##### Footnote: Multiple Plots ########
#######################################
# A Create predictions and 95% CI for 2 df
# A1 Create predictions for 2 degrees of freedom
predBMI <- predict(mod.ns.2, se.fit = TRUE, type = "terms" )
# A2Convert to dataframe 
predBMI <- as.data.frame(predBMI)
# A3 Combine predictions and standard errors
predBMI <- predBMI %>% 
  mutate( pred = fit.ns.avePM.idw..df...2.,
          se = se.fit.ns.avePM.idw..df...2.)
# A4 Compute 95% confidence intervals 
predBMI <- predBMI %>% 
  mutate( lci = pred - 1.96*se,
          uci = pred + 1.96*se)
# A5 Keep only variables we need 
predBMI2 <- predBMI %>% select(pred, se, lci, uci) %>% 
  mutate(Model = "2 df")
# A6 Combine with data 
predBMI.2 <- predBMI2 %>% bind_cols(df)

# B Create predictions and 95% CI for 3 df
# B1 Create predictions for 3 degrees of freedom
predBMI <- predict(mod.ns.3, se.fit = TRUE, type = "terms" )
# B2 convert to dataframe 
predBMI <- as.data.frame(predBMI)
# B3 Combine predictions and standard errors
predBMI <- predBMI %>% 
  mutate( pred = fit.ns.avePM.idw..df...3.,
          se = se.fit.ns.avePM.idw..df...3.)
# B4 Compute 95% confidence intervals 
predBMI <- predBMI %>% 
  mutate( lci = pred - 1.96*se,
          uci = pred + 1.96*se)
# B5 Keep only variables we need
predBMI2 <- predBMI %>% select(pred, se, lci, uci)%>% 
  mutate(Model = "3 df")
# B6 Combine with data 
predBMI.3 <- predBMI2 %>% bind_cols(df)

# C Create predictions and 95% CI for 4 df
# C1 Create predictions for 4 degrees of freedom
predBMI <- predict(mod.ns.4, se.fit = TRUE, type = "terms" )
# C2 convert to dataframe 
predBMI <- as.data.frame(predBMI)
# C3 Combine predictions and standard errors
predBMI <- predBMI %>% 
  mutate( pred = fit.ns.avePM.idw..df...4.,
          se = se.fit.ns.avePM.idw..df...4.)
# C4 Compute 95% confidence intervals 
predBMI <- predBMI %>% 
  mutate( lci = pred - 1.96*se,
          uci = pred + 1.96*se)
# C5 Keep only variables we need
predBMI2 <- predBMI %>% select(pred, se, lci, uci) %>% 
  mutate(Model = "4 df")
# C6 Combine with data 
predBMI.4 <- predBMI2 %>% bind_cols(df)
# C7 combine all the data 
allModels <- bind_rows(predBMI.2, predBMI.3, predBMI.4)

# D Uncenter the data - add mean aveBMI
allModels <- allModels %>% mutate(predBMI = pred + mean(aveBMI),
                                  lciBMI = lci + mean(aveBMI),
                                  uciBMI = uci + mean(aveBMI))

# E Plot 

pdf(paste0(OutputPath, "three_ns_model_plot.pdf"))

# no confidence intervals
ggplot(allModels, aes(avePM.idw)) + 
  geom_line(aes(y = predBMI, color = Model)) + 
  xlab(expression("Average Annual PM"[2.5])) + 
  ylab("Predicted BMI") + 
  ylim(26, 29)

# confidence intervals
ggplot(allModels, aes(avePM.idw)) + 
  geom_line(aes(y = predBMI, color = Model)) + 
  geom_line(aes(y = lciBMI, color = Model), alpha = 0.1) + 
  geom_line(aes(y = uciBMI, color = Model), alpha = 0.1) + 
  xlab(expression("Average Annual PM"[2.5])) + 
  ylab("Predicted BMI (95% CI)") + 
  ylim(26, 29)

dev.off()


