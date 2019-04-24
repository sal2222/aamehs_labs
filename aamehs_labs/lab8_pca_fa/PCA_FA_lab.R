# PCA+FA Lab Session
# Marianthi-Anna Kioumourtzoglou
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# March 13, 2019
# Session 8: PCA+FA
# Special thanks to 
# Yanelli Nunez, Lizzy Gibson, Ahlam Abuawad, and Marianthi Kioumurtzoglou 
# For example code from the 2018 Mixtures Workshop

####***********************
#### Table of Contents ####
# 0:  Preparation 
# 1:  Prepare Exposure Data

# A: Principal Component Analysis 
# 2: Principal Component Analysis 
# 3: Visualize PCA

# B: Exploratory Factor Analysis
# 4: Orthogonal EFA 
# 5: Oblique EFA 
# 6: Compare EFA Model Fit 
# 7: Visualize EFA 
# 8: Assess EFA Scores

####********************
#### 0: Preparation ####
####********************
# 0a Install packages 

# 0a.i Data management packages

#install.packages("tidyverse")
#install.packages("janitor")

# 0a.ii Visualization packages

#install.packages("ggcorrplot")
#install.packages("ggfortify")
#install.packages("gridExtra")
#install.packages("factoextra")
#install.packages("ggrepel")
#install.packages("GPArotation")

# 0a.iii Factor Analysis packages 

#install.packages("psych")

# 0b Load packages

library(tidyverse)
library(dplyr)
library(janitor)
library(lubridate)
library(ggplot2)
library(ggcorrplot)
library(ggfortify)  
library(gridExtra)
library(factoextra)
library(GPArotation)
library(ggrepel)
library(psych)



####*******************************
#### 1: Prepare Exposure Data #####
####*******************************

# 1a Readin Data 

df0 <- read_csv("Boston_PM_constituents.csv")

# 1b Convert date variable to datetime format 

df0$Date <- parse_date_time(df0$Date, "mdy")

# 1c View data 

head(df0)

# 1d Plot potassium concentrations over time

plot(df0$Date, df0$K, type="l")

# 1e Remove days affected by 4th of July 
# fireworks actually emit a significant amount of some constituents, 
# especially potassium 
# so typically air pollution analysts will exclude 4th of July, 
# and the days preceeding it. 

# 1f.i Create a list of days to remove

DaysAffected <- c("7_2", "7_3", "7_4", "7_5")

# 1f.ii Remove those days, for each year 

df0 <- df0 %>% mutate(MonthDay = paste0(month(Date), "_", day(Date))) %>% 
  filter(!(MonthDay %in% DaysAffected )) %>% 
  dplyr::select(-MonthDay)

# 1g Remove Variables
# we  remove PM2.5, since PM2.5 is total mass of particles 
# and represents the sum of mass of the constituent elements
# we will also remove Sodium
# because XRF readings for Na are not so accurate 

df0 <- df0 %>% dplyr::select(-pm25, -Sodium)

# 1h Remove days with missing data 

df0 <- df0 %>% filter(complete.cases(df0))

# 1i Remove non-exposure variables 

df <- df0 %>% dplyr::select(-Date, -ent)

# 1j Check dimensions of dataset

dim(df)

# 1k Summary statistics on dataset

summary(df)

####************************************************
#### A: Principal Components Analysis #####*********
####************************************************

####***************************************
#### 2: Principal Components Analysis #####
####***************************************

# 2a Run PCA 
# scale.: a logical value indicating whether the variables should be scaled
# prcomp requires the variables to have a variance of 1 before doing analyis 
# if scale = TRUE, prcomp will scale our data for us. 

pca.df  <- prcomp(df, scale. = TRUE) 

# 2b Structure of PCA output

ls(pca.df)

# sdev is the scaled singular values
# rotation is the loadings of each component onto each constituent

# 2c Summary of PCA results
# note that the proportion of variance explained is based on the eigenvalues 

summary(pca.df)

# 2c Extract the eigenvalues 
# the eigenvalues are the square of the singular values

eigenvalues.v <- pca.df$sdev^2

# 2d Create table of eigenvalues and percent variance 

# 2d.i Compute percent of variance explained by each principal component

perc_variance.v <- eigenvalues.v/sum(eigenvalues.v)

# 2d.ii Put percent variance in percentage format 

perc_variance.v <- round(100 * perc_variance.v, 1)

# 2d.iii Compute cumulative percent variance explained

cumulative_perc_var.v <- cumsum(perc_variance.v)

# 2d.iv Create dataframe 

eigenvalues.df <- data.frame(1:length(eigenvalues.v), eigenvalues.v, perc_variance.v, cumulative_perc_var.v)

# 2d.v Name dataframe columns 

colnames(eigenvalues.df) <- c("Principal Component", "Eigenvalues", 
                              "Proportion Var Explained", "Cumulative Proportion")

# 2d.vi View eignvalues 

eigenvalues.df

# 2e Assess number of principal components with an eigenvalue > 1
# Since we scaled the exposure data, the variables have a variance of 1 
# so any PC with an eigenvalue > 1 is explaining more of the variability 
# of the data than a column of the original constituent data 

# 7 PC's have eigenvalue >1

###***********************
#### 3: Visualize PCA #####
####***********************

# 3a Scree Plot
# Plots the proportion of variance explained by each component 

# fviz is a specially written function 
# that takes in pca solution objects 
# and then visualizes proportion variance from each PC
# One way to determine the number of factors or components in 
# a data or correlation matrix is to examine
# the “scree" plot of the successive eigenvalues. 
# Sharp breaks in the plot suggest the appropriate number of components
# or factors to extract. 


fviz_eig(pca.df, main = "Percent Variance Explained \n by Principal Component",
         xlab = "Principal component",
         ylim = c(0,40)) 


# 3b Visualization of the loadings 
# Loadings are the weights that each constitient contributes to a component

# 3b.i Extract the rotation or loadings matrix 

loadings.df <- as.data.frame.matrix(pca.df$rotation) 

# 3b.ii Create column with names of constituents

loadings.df$Constituent <- row.names(loadings.df)

# 3b.iii Put loading data in long format 

loadings.long <- loadings.df %>% 
  gather(key = "PC", value = "Loading", -Constituent) 

# 3b.iv Choose just the first 7 principal components 

loadings.long.7PC <- loadings.long %>% 
  filter(PC %in% c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")) 

# 3b.v Plot 

ggplot(loadings.long.7PC, aes(x = Constituent, y = Loading)) + 
  geom_col() +                             # creates a column for each loading
  geom_hline(yintercept = 0, size = 0.2) + # creates a line at 0
  facet_wrap(~ PC) +                       # creates a distinct box for each PC 
  theme_bw() +                             # sets theme options
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = expression("PM"[2.5]~"Constituents"),
       y = "Loadings")

# 3c.i Plot loadings
# Creates a plot showing the loadings for principal component 1 and 2. 

autoplot(pca.df,    # name the pca solution
         data = df, # name the original data
         size = 0.8, colour = 'blue', alpha = 0.5,    
         loadings = FALSE, loadings.colour = 'orange',  
         loadings.label = TRUE, loadings.label.repel = T, 
         loadings.label.size = 2.5, loadings.label.colour = 'black',
         main = "Principal Component Analysis Loading Plot")

# 3c.ii ggplot version 

# ggplot(loadings, aes(x = PC1, y = PC2)) + 
#   geom_point() +
#   geom_label_repel(aes(label = Constituent),
#                    box.padding   = 0.35,
#                    point.padding = 0.5,
#                    segment.color = 'grey50') + 
#   theme_bw() + 
#   theme(legend.position = "bottom") +
#   labs(title = "Variable Loadings on First and Second Factors")

####**********************************************
#### B: Exploratory Factor Analysis #####*********
####**********************************************

# clean the environment 
rm(eigenvalues.df, eigenvalues.v, loadings.df, loadings.long, loadings.long.7PC, perc_variance.v )

# In factor analysis, we need to pre-specify the number of factors in the solution 

# If we have expert knowledge about the variables / exposures 
# we can use that knowledge to inform the number of factors 

# If we do not have expert knowledge, 
# We can use the data to suggest a range of number of factors 
# we can perform a PCA and then either
# use the elbows of the scree plot 
# or the number of components with eigenvalues >1

# In this case, there have been a number of studies of air pollution 
# in Boston, 
# and based on those studies, we expect 5-7 sources. 

####*************************
#### 4: Orthogonal EFA  #####
####*************************

# Orthogonal rotation is used if it is desirable to identify factors 
# that are as independent from one another as possible.
# In Factor Analysis you get values for uniqueness and comparativeness. 

# 4a Create orthogonal EFA solution with 5 factors 

fa_5.ortho <- fa(df,
                 nfactors = 5,            # number of factors in the solution 
                 rotate   = "varimax",    # rotation
                 scores   = "regression", # method to estimate scores
                 fm       = "ml")         # estimation method

# 4b View EFA solution with 5 factors

fa_5.ortho

####***********************
#### 5: Oblique EFA   #####
####***********************

# Oblique rotation (i.e. correlated factors) are commonly used 
# since we often hypothesize our latent variables of interest
# to be correlated with one another.
# Rotated ortogonal solution to get correlated results.
# In the orthogonal solution we only get uncorrelated results. 

# 5a Create oblique solution with 5 factors

fa_5.oblique <- fa(df, 
                 nfactors = 5,             
                 rotate   = "promax",     # promax is a popular oblique rotation
                 scores   = "regression", 
                 fm       = "ml") 

# 5b View oblique solution with 5 factors

fa_5.oblique

# 5c Create oblique solution with 6 factors 

fa_6.oblique <- fa(df, nfactors = 6, rotate = "promax", scores = "regression", fm = "ml") 

# 5d View oblique solution with 6 factors 

fa_6.oblique

# 5e Create oblique EFA solutions with 7 factors  

fa_7.oblique <- fa(df, nfactors = 7, rotate = "promax", scores = "regression", fm = "ml") 

#####*****************************
#### 6: Compare EFA Model Fit  #### 
#####*****************************
# If we do not have expert knowledge to guide us in choosing the best factor solution 
# we can used eBIC to compare how well models fit the data 

# from the documentation"
# "eBIC -- When normal theory fails (e.g., in the case of non-positive definite matrices), 
# it useful to examine the empirically derived eBIC
# based upon the empirical chi^2 - 2 df."

# 6a Create table of EFA models and their eBIC scores

# 6a.i Construct Table

fit.efa <- as.data.frame(rbind(cbind("5 Factor", "Promax", round(fa_5.oblique$EBIC)),
                               cbind("6 Factor", "Promax", round(fa_6.oblique$EBIC)),
                               cbind("7 Factor", "Promax", round(fa_7.oblique$EBIC))))

# 6a.ii Rename columns

names(fit.efa) <- c("Model", "Rotation", "EBIC")

# 6b Examine fit table 

fit.efa

# Choose Promax 7 factor model based on fit statistics and interpretability.

#####*********************
#### 7: Visualize EFA  #### 
#####*********************

# 7a Extract loadings

# 7a.i Create table of loadings 

loadings <- data.frame(fa_7.oblique$loadings[])

# 7a.ii Create column with names of constituents

loadings$Constituent <- rownames(fa_7.oblique$loadings)

# 7c Identify the factor with the largest contribution to each constituent

loadings$Max <- colnames(loadings[,1:7])[max.col(loadings[,1:7], ties.method = "first")]

# 7d Plot Loadings 

# 7d.i Rename factor columns 
loadings <- loadings %>% 
  rename(Factor1 = ML1, 
         Factor2 = ML2, 
         Factor3 = ML3, 
         Factor4 = ML4, 
         Factor5 = ML5, 
         Factor6 = ML6, 
         Factor7 = ML7)

# 7d.i Put loading data in long format 

loadings.long <- loadings %>% 
  dplyr::select(-Max) %>% 
  gather(key = "Factor", value = "Loading", -Constituent)

# 7d.ii Choose just the first 7 factors
# here we will plot all 7 factors, but you can use this line 
# to choose to plot fewer factors

loadings.long.6F <- loadings.long %>% 
  filter(Factor %in% c("Factor1", "Factor2", "FactorL", "Factor4", "Factor5", "Factor6")) 

# 7d.iii Loadings plot 

ggplot(loadings.long.6F, aes(x = Constituent, y = Loading)) + 
geom_col() +
geom_hline(yintercept = 0, size = 0.2) +
facet_wrap(~ Factor) + 
theme_bw() + 
theme(strip.background = element_rect(fill = "white")) +
labs(x = expression("PM"[2.5]~"Constituents"),
     y = "Loadings")


# 7e Plot Loadings
ggplot(loadings, aes(x = Factor1, y = Factor2)) + 
geom_point() +
geom_label_repel(aes(label = Constituent),
                   box.padding   = 0.35,
                   point.padding = 0.5,
                   segment.color = 'grey50') + 
theme_bw() + 
theme(legend.position = "bottom") +
labs(title = "Variable Loadings on First and Second Factors")

#####**************************
#### 8: Assess EFA Scores  #### 
#####**************************

# Each day has a score for each factor 
# The score is the degree to which that factor contributed to 
# that day's PM

# 8a Extract Scores

# 8a.i Create table of scores

scores <- data.frame(fa_7.oblique$scores[])

# 8a.ii Create column with names of constituents

scores$Constituent <- rownames(fa_7.oblique$scores)

# 8b For each day, identify the factor that contributed the most to the PM2.5

scores$Max <- colnames(scores)[max.col(scores, ties.method = "first")]

# 8c Add day column 

scores$Date <- df0$Date

# 8c View scores table 
head(scores)

