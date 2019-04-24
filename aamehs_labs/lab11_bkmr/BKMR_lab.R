# BKMR Session
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# April 10, 2018
# Session 10: Bayesian Kernel Machine Regression
# Thanks to Jennifer F. Bobb at Kaiser Permanente
# for creating the BKMR package 
# and her online tutorial at 
# https://jenfb.github.io/bkmr/overview.html

##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Load Data 
# 2: Fit BKMR Model
# 3: BKMR Convergence
# 4: Plot BKMR Results

####********************
#### 0: Preparation ####
####********************
# 0a Install packages 

install.packages("bkmr")

# 0b Load packages

library(readr)
library(dplyr) 
library(bkmr)
library(ggplot2)

#0c Declare directories


###############################################################
####*******************
##### 1: Load Data ####
####*******************
# 1a Load data 

df <- read_csv("county_bmi_gases_confounders_10.csv")

# 1b Review data

head(df)

summary(df)

View(df)

####***********************
#### 2: Fit BKMR Model ####
####***********************

# 2a Set Seed

set.seed(111)

# 2b Create matrices of the expsoures and confounders
# bkmr requires matrices for the exposure matrix

gases <- cbind(df$aveo3.idw, df$aveco.idw, df$aveso2.idw)

confounders <- as.matrix(df[,6:15])

# Please note that we are not including here Region which is a categorical variable
# BKMR requires continuous confounders
# We can transform the Regions into separate 8 dummy variables and include 7 of the
# 8 dummies in the model


# 2c Fit BKMR model

fitkm <- kmbayes(y = df$aveBMI,   # outcome variable
                 Z = gases,       # matrix of exposures
                 X = confounders, # matrix of confounders
                 iter = 10000,    # number of iterations
                 verbose = FALSE, # print interim progress of model?
                 varsel = TRUE) 

# 2d Save BMKR Model

saveRDS(fitkm, "initial_bkmr_solution.RDS")

# fitkm <- readRDS("initial_bkmr_solution.RDS")

####**************************
#### 3: BKMR Convergence  ####
####**************************
# BKMR estimates associations by running a series of iterations 
# (think conceptually similar to k-means). 
# the model achieves convergence when the estimates are stable. 
# when we implement BKMR, we check for convergence of the model 
# If we do not have convergence, then we can increase the number of iterations 
# or take other steps to improve model convergence 

# 3a Review trace plots 
# one way we can review convergence is by looking at the trace plots 
# trace plots show the behavior of the model at each iteration for the estimated parameters

# a trace plot with flat sections or steep jumps
# indicates that the model did not reach convergence 
# when these are towards the later iterations. 
# we want to see stability at the ends

# (also note that we asked for 10K iterations and we only plot 5K -- "burn in")

# the betas are the coefficient for each confounder

TracePlot(fit = fitkm, par = "beta")
TracePlot(fit = fitkm, par = "beta", comp = 1)

TracePlot(fit = fitkm, par = "beta", comp = 2)
TracePlot(fit = fitkm, par = "beta", comp = 10)

# "sigsq.eps" and "r" are tuning parameters of the bayesian model 
# for the purposes of this class we cannot go into details
# for those interested in more, please feel free to check
# Jen Bobb's tutorial: https://jenfb.github.io/bkmr/overview.html
# and her paper, and let us know if you have any questions!

TracePlot(fit = fitkm, par = "sigsq.eps")

TracePlot(fit = fitkm, par = "r", comp = 1)
TracePlot(fit = fitkm, par = "r", comp = 2)


# 3b Standardize confounders to help convergence

confounders2 <- scale(df[,6:15])

# 3c Run model with standardized covariates 

fitkm2 <- kmbayes(y = df$aveBMI,   
                  Z = gases,        
                  X = confounders2, 
                  iter = 10000,     
                  verbose = FALSE,  
                  varsel = TRUE) 

# fitkm2 <- readRDS(paste0(DataPath, "second_bkmr_solution.RDS"))

# 3d Review trace plots

TracePlot(fit = fitkm2, par = "beta")
TracePlot(fit = fitkm2, par = "beta", comp = 2)
TracePlot(fit = fitkm2, par = "beta", comp = 10)

TracePlot(fit = fitkm2, par = "sigsq.eps")

TracePlot(fit = fitkm2, par = "r", comp = 1)
TracePlot(fit = fitkm2, par = "r", comp = 2)

# 3e Save BMKR Model

saveRDS(fitkm2, "second_bkmr_solution.RDS")

####**************************
#### 4: Plot BKMR Results ####
####**************************

# 4a Readin BKMR model 

fitkm <- readRDS("second_bkmr_solution.RDS")

# 4b Posterior inclustion probabilities 
# The PIP is the probability that the exposure is an important 
# part of the mixture

ExtractPIPs(fitkm, z.names = c("ozone", "carbon monoxide", "sulfur dioxide"))

# 4b Plot Univariate exposure-response

# 4b.i Estimate exposure-reponse for each exposure 
# When the other exposures are at their median value 
# In the presence of interaction, 
# there is not a single exposure-response curve
# the association between the exposure and outcome depends on the other exposures

pred.resp.univar <- PredictorResponseUnivar(fit = fitkm, 
                                            q.fixed = 0.5, # set other exposures to 50th percentile
                                            z.names = c("ozone", "carbon monoxide", "sulfur dioxide")) 

# 4b.ii Create Plot 

ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = "identity") + 
  facet_wrap(~ variable, scales = "free_x") +
  ylab("Change in Expected aveBMI") + 
  xlab("Gas Concentration")


# 4c Plot exposure-response curves at various levels of co-exposures

# 4c.i Estimate exposure-response for a range of co-exposure levels

pred.resp.bivar <- PredictorResponseBivar(fit = fitkm, min.plot.dist = 1)

# 4c.ii Extract exposure-response and their uncertainity
# for specific levels of co-exposure

pred.resp.bivar.levels <- PredictorResponseBivarLevels(
  pred.resp.df = pred.resp.bivar, # the predictions from which to extract
  Z = gases,                      # matrix of exposures
  qs = c(0.1, 0.5, 0.9))          # the quantiles of the co-exposures

# 4c.iii Reorganize exposure-response values 

pred.resp.bivar.levels <- pred.resp.bivar.levels %>% 
  mutate(variable1 = case_when(
    variable1 == "z1" ~"ozone", 
    variable1 == "z2" ~"carbon monoxide", 
    variable1 == "z3" ~"sulfur dioxide"  )) %>% 
  mutate(variable2 = case_when(
    variable2 == "z1" ~"ozone", 
    variable2 == "z2" ~"carbon monoxide", 
    variable2 == "z3" ~"sulfur dioxide"  ))

# 4c.iv Create Plot

ggplot(pred.resp.bivar.levels, aes(z1, est)) + 
  geom_smooth(aes(col = quantile), stat = "identity") + 
  facet_grid(variable2 ~ variable1, scales = "free_x") +
  ggtitle("Exposure-Response (Exposure 1 | quantiles of Exposure 2)") +
  xlab("Concentration of Exposure 1")


# the last row looks only blue b/c the lines are completely overlapping!

# 4d Plot exposure-response of total mixture

# 4d.i Estimate exposure-response for each quantile 
# estimates exposure-response if each exposure is set to that quantile 
# ie the first value is the response when 
# ozone is at its tenth quantile of concentration (0.04)
# co is at its tenth quantile of concentration  (8.3)
# so2 is at its tenth quantile of concentration (0.5)

# The estimates are the change in expected aveBMI 
# as the total mixture changes from 
# the concentrations at the q.fixed level 
# (in this example, the median concentration)
# to that quantile of exposure

responses.overall <- OverallRiskSummaries(fit = fitkm,     # model 
                                      y = df$aveBMI,   # outcome 
                                      Z = gases,       # exposures
                                      X = confounders2, # matrix of confounders 
                                      qs = seq(0.25, 0.75, by = 0.05), # vector of quantiles
                                      q.fixed = 0.5, 
                                      method = "exact") # estimation method

# 4d.ii Table of Responses

responses.overall 

# 4d.iii Plot

ggplot(responses.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange()


# 4d.iv Set q.fixed to 0.25
# This is just to illustrate the role of the q.fixed arguement

responses.overall <- OverallRiskSummaries(fit = fitkm,      # model 
                                          y = df$aveBMI,    # outcome 
                                          Z = gases,        # exposures
                                          X = confounders2, # matrix of confounders 
                                          qs = seq(0.25, 0.75, by = 0.05), # vector of quantiles
                                          q.fixed = 0.25, 
                                          method = "exact") # estimation method

# 4d.iii Plot

ggplot(responses.overall, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange()


