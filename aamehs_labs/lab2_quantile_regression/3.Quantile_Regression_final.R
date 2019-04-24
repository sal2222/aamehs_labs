# Quantile Regression Session
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# Dec 20, 2018
# Session 2: Quantile Regression 
# Step 3: Session
##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Load Data 
# 2: Median Model 
# 3: Compare Quantiles 
# 4: Group Exercises 
# Footnote: Standard Errors
# Footnote: ggplot2 options
########################
#### 0: Preparation ####
########################
# 0a Install packages 

# quantreg features the quantile regression and associated functions
# install.packages("quantreg")

# ggplot2 has plotting commands 
# install.packages("ggplot2")

# 0b Load packages
library(tidyverse) 
library(quantreg)


#0c Declare directories

DataPath   <- "/Users/slewa/Desktop/AAMEHS/lab2_quantile_regression/"
OutputPath <- "~/Desktop/AAMEHS/lab2_quantile_regression/outputs/" 
# 
# DataPath   <- paste0("/Users/marianthi_anna/Dropbox/AAMEHS_Spring2019/Labs/Session2_QuantileRegression/")
# OutputPath <- paste0("/Users/marianthi_anna/Dropbox/AAMEHS_Spring2019/Labs/Session2_QuantileRegression/")
# 
# 
###############################################################
#######################
##### 1: Load Data ####
#######################
# 1a load data 

df <- read_csv(paste0(DataPath, "county_bmi_pm_confounders_10.csv"))

# 1b remove the variables we will not need 

df <- df %>% select(-AREA, -PERIMETER, -geometry)

##########################
##### 2: Median Model ####
##########################

# 2a Create median model 
# how do we tell the rq function which quantile to estimate? 

?rq 

# the parameter tau is the quantile(s) to be estimated 
# median is 50th percentile, so tau = 0.5

mod50 <- rq(aveBMI ~ avePM.idw + maleunemp  + femaleunemp + ltHS + medhinc +
              medhval + per.black + per.latinx + per.asnam + climate_region, 
            data = df, 
            tau = 0.5)

# Class Question ######################################
# Do we have any warnings? What do they mean? How do we learn more? 
####################

FAQ()

# 2b Median model results
# we can use summary just as before 
# alpha sets the Type I error rate for the confidence intervals

summary(mod50, alpha = 0.05)

# However, the output is slightly different 
# We get upper and lower bounds rather than standard errors and T tests. 
# I explain how to get standard errors in the Footnotes Section

# Class Question ######################################
# How do we interpret the coefficient for avePM.idw? 
#**********************

# 2c Compare to model for mean aveBMI

modMean <- lm(aveBMI ~ avePM.idw + maleunemp  + femaleunemp + ltHS + medhinc +
                medhval + per.black + per.latinx + per.asnam + climate_region, 
              data = df)

# Let's compare the coefficients of the models 

mod50.coeff   <- summary(mod50, alpha = 0.05)$coefficients[,1]
modMean.coeff <- summary(modMean)$coefficients[,1]
coeff.table0  <- data.frame(mod50.coeff, modMean.coeff)

# Let's make the table more readable
coeff.table <- coeff.table0 %>% 
  #compute percentage difference
  mutate(percent.diff = 100*(mod50.coeff - modMean.coeff)/modMean.coeff) %>% 
  # round numbers 
  mutate(mod50.coeff   = round(mod50.coeff, 3), 
         modMean.coeff = round(modMean.coeff, 3),
         percent.diff   = round(percent.diff, 1))

# Class Question ######################################
# Are the models similar? 
#**********************


# 2d Compare models with plots 
# for these plots, we will use ggplot2 


# 2d.i combine estimates from each model into a single dataframe 
# assemble estimates from each model
# note that the code here is different for the rq model and the lm model 
# since the summary.rq() directly outputs confidence intervals 
# and summary.lm outputs standard errors.

# FYI, the : notation calls for each element between 1 and 3.
# including 1 and 3
# and c() makes a vector out of those three elements

MedianModel <- c(summary(mod50, alpha = 0.05)$coefficients[2,1:3])

coeff.lm    <- summary(modMean)$coefficients
MeanModel <- c(coeff.lm[2,1], 
               coeff.lm[2,1] - 1.96 * coeff.lm[2,2], 
               coeff.lm[2,1] + 1.96 * coeff.lm[2,2])

# create dataframe 

coeff.table <- rbind(MedianModel, MeanModel)
coeff.table <- as.data.frame(coeff.table, stringsAsFactors = FALSE)

# set names for dataframe

names(coeff.table) <- c("coeff", "lci", "uci")
coeff.table        <- coeff.table %>% 
                      mutate(ModelName = c("Median Model", "Mean Model"))

# 2d.ii Forest plot
# with ggplot, we can keep our plot in the environment
# as a gpglot object 
# and then display it in the plot panel 
# by entering the plot's name

             # defines what dataset ggplot will use
fp.mean.median <- ggplot(data = coeff.table,     
             # aes() defines which variables the geoms will use   
             aes( # defines variable for the x axis
                 x = ModelName,  
                 # defines the variable for the point along the y axis
                 y = coeff,      
                 # defines the lower bound of the confidence interval
                 ymin = lci,     
                 # define the upper bound of the confidence interval 
                 ymax = uci)) +  
  # creates a point (y) with line defined by ymin and ymax
  geom_pointrange() +   
  # creates lines with bars, i.e. here the CIs
  geom_errorbar()+      
  # add a dashed line at y=0
  geom_hline(aes(yintercept = 0.0), lty = 2) +
  # labels for axes
  xlab("Model Name") +    
  ylab(expression("Coefficient for PM"[2.5]~" (95% CI)"))

fp.mean.median # produces the plot in the plots panel

###############################
##### 3: Compare Quantiles ####
###############################
# often, we are interested in how associations vary for different quantiles 

# 3a Two tau's 
# the quantreg package incldues a way to estimate 
# models for multiple tau's at once 
# we create a vector of the tau's

mods25.50 <- rq(aveBMI ~ avePM.idw + maleunemp  + femaleunemp + ltHS + medhinc + 
                  medhval + per.black + per.latinx + per.asnam + climate_region, 
                data = df, tau= c(0.25, 0.50)) # c() creates a vector

summary(mods25.50)

# 3b Plot two quantile models 
# we will follow the same procedure as before. 

# 3b.i Combine estimates from each model into a single dataframe 
# assemble estimates from each model
# we will save the summary as an object 
# so we do not have to keep re-summarizing the models. 

summary25.50 <- summary(mods25.50, alpha = 0.05)

# we use the brackets to specify which model we are extracting

Model25th   <- c(summary25.50[[1]]$coefficients[2,1:3])
Model50th   <- c(summary25.50[[2]]$coefficients[2,1:3])

# create dataframe 

coeff.table <- rbind(Model25th, Model50th)
coeff.table <- as.data.frame(coeff.table, stringsAsFactors = FALSE)

# set names for dataframe

names(coeff.table) <- c("coeff", "lci", "uci")
coeff.table        <- coeff.table %>% 
                      mutate(ModelName = c("Model 25th", "Model 50th")) 

# 3b.ii Plot

fp25.50 <- ggplot(data=coeff.table, # defines what dataset we are using
             aes(x=ModelName,  # defines variable for the x axis
                 y=coeff,      # defines the variable for the point along the y axis
                 ymin=lci,     # defines the lower bound of the confidence interval
                 ymax=uci)) +  # define the upper bound of the confidence interval   
  geom_pointrange() +          # creates a point (y) with line defined by ymin and ymax        
  geom_errorbar()+             # creates lines with bars
  geom_hline(aes(yintercept=0.0), lty=2) + # add a dashed line at y=0 
  xlab("Model Name") +         # labels for axes
  ylab(expression("Coefficient for PM"[2.5]~" (95% CI)"))

fp25.50

# 3c Plot Multiple Quantiles
# we can just use the same procedure as before, 
# but with additional models for every quantile of interest

# 3c.i Create the models
# seq() creates a sequence with intervals set by the by arguement 

TauList <- seq(0.1, 0.9, by = 0.1)
TauList

qr.mods  <- rq(aveBMI ~ avePM.idw + maleunemp  + femaleunemp + ltHS + medhinc + 
               medhval + per.black + per.latinx + per.asnam + climate_region, 
               data = df, 
               tau = TauList)

# 3c.ii Assemble estimates from each model

summary.qr.mods <- summary(qr.mods, alpha = 0.05)

Model10th   <- c(summary.qr.mods[[1]]$coefficients[2,1:3])
Model20th   <- c(summary.qr.mods[[2]]$coefficients[2,1:3])
Model30th   <- c(summary.qr.mods[[3]]$coefficients[2,1:3])
Model40th   <- c(summary.qr.mods[[4]]$coefficients[2,1:3])
Model50th   <- c(summary.qr.mods[[5]]$coefficients[2,1:3])
Model60th   <- c(summary.qr.mods[[6]]$coefficients[2,1:3])
Model70th   <- c(summary.qr.mods[[7]]$coefficients[2,1:3])
Model80th   <- c(summary.qr.mods[[8]]$coefficients[2,1:3])
Model90th   <- c(summary.qr.mods[[9]]$coefficients[2,1:3])

# create dataframe 

coeff.table <- rbind(Model10th, Model20th, Model30th, Model40th, 
                     Model50th, Model60th, Model70th, Model80th, 
                     Model90th)

coeff.table <- as.data.frame(coeff.table, stringsAsFactors = FALSE)

# set names for dataframe

names(coeff.table) <- c("coeff", "lci", "uci")
coeff.table        <- coeff.table %>% 
          mutate(ModelName = c("Model 10th", "Model 20th", "Model 30th", "Model 40th",
                               "Model 50th", "Model 60th", "Model 70th", "Model 80th", 
                               "Model 90th"))

# 3b.ii Plot

fp.qr.mods <- ggplot(data=coeff.table, # defines what dataset we are using
                  aes(x=ModelName,  # defines variable for the x axis
                      y=coeff,      # defines the variable for the point along the y axis
                      ymin=lci,     # defines the lower bound of the confidence interval
                      ymax=uci)) +  # define the upper bound of the confidence interval   
  geom_pointrange() +               # creates a point (y) with line defined by ymin and ymax        
  geom_errorbar()+                  # creates lines with bars
  geom_hline(aes(yintercept=0.0), lty=2) + # add a dashed line at y=0 
  xlab("Model Name") +              # labels for axes
  ylab(expression("Coefficient for PM"[2.5]~" (95% CI)"))

fp.qr.mods



#############################
#### 4: Group Exercises  ####
#############################

# 4a Estimate the adjusted association between county-average annual PM2.5 
# and the 35th percentile of county-average BMI 

# 4b Plot the coefficient and confidence intervals for avePM.idw for the models 
# for the 35th and 65th percentile

# 4c Plot the coefficient and confidence interval for the association between 
# medhval and 75th percentile aveBMI


####################################
#### Footnote: Standard Errors  ####
####################################

# The author of the rq package, Roger Koenker, 
# included multiple possible algorithms for estimating 
# the uncertainity of the model. 
# more information can be found with 
# ?summary.rq
# under the heading 'se'
# Briefly, when the data has fewer than 1001 observations, then 
# the model will use the rank method, which does not estimate standard errors, 
# only confidence intervals. 
# For larger datasets, or if you set the se manually, 
# the model will use a different method that does compute standard errors. 
# In that case, you can use similar code as with linear regression 
# to extract confidence intervals. 



######################################
#### Footnote 2: ggplot2 options  ####
######################################

# More ggplot2 options
# there are many parts of the plot you can specify 
# if you don't like the defaults 
# check the ggplot cheatsheet or just look online 
# for help 
# here is just some example code. 

fp.fancy <- ggplot(data=dgk, 
             # tell geom_functions which variables to use
             aes(x=ModelName, y=coeff, ymin=lci, ymax=uci)) + 
  # point is a diamond , increase size 
  geom_pointrange(shape = 18, size = 1.5) +       
  # increase the size of the error bars
  geom_errorbar(size = 1.5)+ 
  # changes the color of the line
  geom_hline(aes(yintercept=0.0), lty=2, color ="grey") + 
  # flip coordinates (puts labels on y axis)
  coord_flip() +                                     
  xlab("Model\nName") + ylab(expression("Coefficient for PM"[2.5]~" (95% CI)")) +
  # change the angle of the title of the y-axis
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +  
  # change the size of the axis titles
  theme(axis.title = element_text(size = 28)) +                 
  # change the size of the axis text
  theme(axis.text = element_text(size = 20)) +      
  # use a white background without gridlines
  theme(panel.background = element_rect(fill = 'white', color = "black")) 

fp.fancy

