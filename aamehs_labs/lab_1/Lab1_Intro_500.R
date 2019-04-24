# Intro to R
# Session 1
# Advanced Analytic Methods for Environmental Epidemiology
# Jan 23, 2019
##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Examine Objects
# 2: Data Manipulation 
# 3: Summarizing Data 
# 4: Group Exercises
# 5: Plotting Data
# 6: Linear Regression
# 7: Group Exercises
##################################################
########################
#### 0: Preparation ####
########################

# 0a load packages
library(tidyverse)


# 0b declare folder paths 
# check the intro_to_R.docx document if you have difficulty finding the file path for your computer

DataPath   <- "/Users/slewa/Desktop/AAMEHS/lab_1/"
OutputPath <- "~/Desktop/AAMEHS/lab_1/outputs/" # Outputs includes plots, summary stats, models


# 0b read in data 

df       <- read_csv(paste0(DataPath, "county_pm_bmi_census_500.csv"))

# read_csv() is a bit faster than read.csv, and has more predictable defaults. 
# even RStudio uses read_csv as the default
# read_csv() comes from the readr package

##############################
#### 1: Examine Objects #####
##############################

# ls(), aka list, shows you  what objects you have saved in R

ls()

# Examine the R data.frame.  
# The df data.frame is an object whose dimensions are 972 observations by 14 variables.  
# Typically each row is an observation and each column is a variable. 

dim(df) # size of the object 
names(df) # names of columns
head(df) # show first 6 rows 
View(df) # shows data as cells


#############
# Variable Codebook 
# FIPS - identifier for each county 
# aveBMI- average BMI for each county in 2010, based on the BRFSS Survey from US CDC 
# avePM.idw - annual ambient PM2.5, based on inverse-distance weighted average of EPA monitors
# climate.region - from NOAA, divides county in area with similar climates
# 2010 American Community Survey of US Census
# per.black, per.latinx, per.asam, per.white - racial composition 
# med.hinc - median household income 
# med.hval - median household vlaue among homeowners 
# lt.hs - percentage of population over 18 with less than high school education 
# female.unemp, male.unemp - unemployment rate by sex
# unemployment has a specific definition, not just 'not working', must be seeking work. 

#############

# Lets examine a single column in the dataframe. 
# Columns in a data.frame may be directly accessed using 
# the data.frame$variable method.
# check the base R cheatsheet and chapter 3 of Paradis tutorial to review what mode class and length mean 

mode(df$avePM.idw) # can be numeric, integer, character, logical
length(df$avePM.idw)

# in R, we can include . and _ in our variable names 
# this is not always the case in other programs: for exmaple, in SAS the . has other uses

# if you do not specify the dataframe, R will assume that avePM.idw is its own object in the environment 
# and will not be able to find it

mode(avePM.idw)

# dataframes are also a grid, so we can look/summon/extract values based on their row+column position 
# value <- dataframe[row, column]

df[1,2]

# we can also extract rows or columns this way 

row1    <- df[1,]
column1 <- df[,1]

rm(row1, column1)
# c() is a function that allows you to combine values into a vector 

a <- c(1,2,3)
a
a + 5

###############################
#### 2: Data Manipulation #####
###############################

# often our raw data aren't ready for analysis, and we need to remove observations, create variables, 
# combine datasets, etc.

# 2a Creating Columns 
# we can create new columns by assigning them 

df$aveBMIsq <- df$aveBMI * df$aveBMI 


df1 <- mutate(df, aveBMIsq = aveBMI*aveBMI)

# within mutate(), the function knows that you are referring to columns withing df, so you do not 
# need the df$ notation 
# we can also use if else statements within dplyr 
# dplyr automatically does operations row by row 
# so adding two columns a and b would create a third column c with values 
# c1 = a1+ b1, c2=a2+b2, c3 = a3+b3

# We often create new dataframes when manipulating data 
# It is important to use consistent naming so you can keep track of what you are doing. 
# If we use the same name for the dataframe, then we are rewriting that dataframe, and we lose previous versions 
# Often, we want to use sequential of informative names so that we can keep track of our progress 
# and catch errors 
# sequential names: df1, df2, df3
# informative names: df_onlyobese, df_bmisq
# naming might seem trivial, but it has helped me stay organized and fix my mistakes
# normally I don't create new dataframe for every transformation, 
# but in this code I did so that is it easy to follow how the code changes the data

# if_else()
# if_else(logical statement, value if true, value if false)
df2 <- mutate(df1, aveObese = if_else(aveBMI >= 30, "average obese", "average not obese"))

# the cut() function 
# the "cut()" function will divide a variable into a set number of levels. 
# We will divide aveBMI into five levels, and then store that information in a new variable, aveBMI5
# The basic cut() function chooses five equally spaced intervals of aveBMI. 
# Since cut() creates intervals that are open on the lower side and closed on the other 
# In the lowest interval, the lowest value will not be included
# since it is the "open end point" 
# the arguement include.lowest = TRUE will include the lowest value in the lowest bin 
# we almost always want to inlcude the lowest value. 

df3 <- mutate(df2, aveBMI5 = cut(aveBMI, 5), include.lowest=TRUE)

# another way to cut your data is to specify cut points that are quintiles of aveBMI. 
# The categories have equal numbers of subjects, but are not equally spaced on the aveBMI axis.  

df2$aveBMI5 <- cut(df2$aveBMI, quantile(df$aveBMI, c(0,.2,.4,.6,.8,1)), include.lowest=TRUE)
df3         <- mutate(df2, aveBMI5 =  cut(aveBMI, quantile(aveBMI, c(0,.2,.4,.6,.8,1)), include.lowest=TRUE))

# Note that here in line 152 I create a new dataframe, df3, and rewrite in in line 158 
# These lines are making the same transformation so its not an issue to rewrite
# In the first draft of code, I usually avoid rewriting objects 
# until I am certain that that line works correctly

# dplyr provides a special syntax to pass along objects

df4 <- df3 %>% mutate(aveBMIsq = aveBMI*aveBMI)

# same result as
df4.a <- mutate(df3, aveBMIsq = aveBMI*aveBMI)

# these pipes can be connected, like this 

df4 <- df3 %>% mutate(aveBMIsq = aveBMI*aveBMI) %>% 
               mutate(aveObese = if_else(aveBMI >=30, "average obese", "average not obese"))

# 2b Subsetting Data 
# we can use the filter() command to keep only rows/observations we are interested in 

df.obese <- df4 %>% filter(aveObese == "average obese")

# filter() can accept a range of logical statements 

df.obese2 <- df4 %>% filter(aveBMI >=30)

# similarly, we can use the select() command to keep only the columns/variables we are interested in 

df.obese3 <- df4 %>% select(fips, aveObese)

# we can also use select() to remove columns with the minus sign 
# here we are removing these three spatial variables which helped us connect the pm and county data, but are no longer necesary 

df5 <- df4 %>% select(-AREA, -PERIMETER, -geometry)


# 2c Combining Data 
# we can combine data with the join commands 
# the dplyr cheatsheet is a good reference for understanding the different ways of combining data. 
# here we will use left_join() to add unemployment data from the 2010 US census 

unemp <- read.csv(paste0(DataPath, "2010_census_unemployment_data.csv"))

# we join by fips codes, so we combine rows that have matching fips

df6 <- left_join(df5, unemp, by = "fips")

# we run into an issue: 
# "incompatible types" 
# which column is an integer, and which one is character? 

is.integer(df5$fips)
is.character(df5$fips)
is.integer(unemp$fips)
is.character(unemp$fips)

# let's convert the fips column in unemp to integer 

unemp.fipschara <- unemp %>% mutate(fips = as.character(fips))

# now let's try the join again 

df6 <- df5 %>% left_join( unemp.fipschara, by = "fips")

# look at the new dataframe... does it contain the unemployment variables? 

View(df6)

# There are other join options 
# for example, inner_join, will only keep observations with a successful match, 
# and full_join will keep all observations from both data frames 
# check the dplyr cheatsheet to see your options
# CRAN hosts several excellent cheatsheets on its site 
# https://www.rstudio.com/resources/cheatsheets/
# I still refer to my cheatsheets

# for example, inner_join() only keeps rows that do have a match for fips. 

df.inner <- df5 %>% inner_join( unemp.fipschara, by = "fips")

# Here is another example of a join that is part of a pipe 
# I will cover this sort of code in the data manipulation workshop 

v   <- read_csv(paste0(DataPath, "noaa_climate_regions_fips.csv"))
v   <- v %>% mutate(state_fips_code = as.character(state_fips_code))

df7 <- df6 %>% 
  mutate(state_fips_code = stringr::str_sub(fips,0,2)) %>% 
  left_join(v, by = "state_fips_code") %>% 
  select(-state_fips_code)

# 2d Cleaning the environment 
# let's remove dataframe we are no longer using
rm(column1, df, df.inner, df.obese, df.obese2, df.obese3, 
   df1, df2, df3, df4, df5, df6, row1, unemp, unemp.fipschara, v, a)
###############################
#### 3: Summarizing Data  #####
###############################

# We ALWAYS need to look at our data before starting our analyses 
# summarizing data is an easy way to catch errors 
# E.g., the first time I created this data set, the maximum average BMI was 400,
# and I realized my scale was wrong 

# just as with data manipulation, there is a base R and a tidyverse method for many tasks 
# Personally, I mostly use the tidyverse method. 
# the tidyverse method outputs dataframes, which I can then easily save or use in other analyses 
# the summarize command also allows me to compute multiple summary stats, for multiple variables, at the same time. 
# the base R commands are simple and require less typing 
# I will present both here. 

# 3a Means
# we can compute the mean of any vector 

MeanPM1 <- mean(df7$avePM.idw)
MeanPM2 <- df7 %>% summarize(meanPM = mean(avePM.idw))

# what is the different between MeanPM1 and MeanPM2? 

class(MeanPM1)
class(MeanPM2)

# 3b Medians 
# we can compute the median of any vector 

MedianPM1 <- median(df7$avePM.idw)
MedianPM2 <- df7 %>% summarize(medianPM = median(avePM.idw))

# 3c Quantiles 

QuantilePM1 <- quantile(df7$avePM.idw, c(0.1,0.9))
QuantilePM2 <- df7 %>% summarize(tenthPM = quantile(avePM.idw, 0.1),
                                 ninetithPM = quantile(avePM.idw, 0.9))

# do two methods yield the same results? 

identical(QuantilePM1[1], QuantilePM2$tenthPM[1])
QuantilePM1[1] 
QuantilePM2$tenthPM[1]


# the dplyr cheat sheet shows the list of the many functions you can use with summarize 

# 3d Multiple summaries 
# in base R, we can quickly create summaries with the summary and table commands 
# you can get the summary of an individual variable 

summary(df7$avePM.idw)
summary(df7$aveBMI)

# the results for summary are not very interesting for character variables 

summary(df7$aveObese) 
table(df7$aveObese) # table is more useful for character variables

# we can rewrite the aveObese variable as a factor
# the results are more useful for factors

df8 <- df7 %>% mutate(aveObese = as.factor(aveObese))
summary(df8$aveObese)

# you can get the summary of every variable in a dataframe 

summary(df8) 

# you can create count cross-tabulations with the table() command 
# aka frequency tables

table(df8$aveObese, df8$state)

# we can save the results 

ObesityByState1 <- table(df8$aveObese, df8$state)

# One advantage of these commands is that they are simple - 
# you don't have to write as much as with the dplyr method 
# the main disadvantage is that the result is not a dataframe, and can be harder to manage. 

# Multiple summary stats in dplyr 
# we can create cross-tabulations using the group_by function and n()
# n() here does not require its own arguement

df8 %>% group_by(state, aveObese) %>%
  summarize(Count = n())

# we can save the results 

ObesityByState2 <- df8 %>% group_by(state, aveObese) %>%
  summarize(Count = n())

# What is the highest county-average PM in each state? 
# (among the counties measured!) 

PM25ByState <- df8 %>% group_by(state) %>% 
  summarize(MaxPM = max(avePM.idw))

#3e  Pipes 
# using the %>% pipes, we can combine manipulations into a single step

# what is the average PM and BMI within Pennsylvania? 

MeanPM_MeanBMI <- df8 %>% filter(state == "PA") %>% 
  summarize(meanPM = mean(avePM.idw), 
            meanBMI = mean(aveBMI))

#3f Correlations 
# we can compute the pearson or spearman correlation with cor() 

cor(df8$aveBMI, df8$avePM.idw, method = "pearson")
cor(df8$aveBMI, df8$avePM.idw, method = "spearman")

###############################
#### 4: Exercises         #####
###############################

# Data Manipulation 
# 4a create a categorical variable for annual average PM2.5 with 4 quantiles 
# 4b create dummy variable for whether a county is in New York, such that 1= in NY and 0 = not NY 
# 4c create dummy variable for whether a county's annual average PM2.5 is greater or smaller than the mean PM2.5 across all counties
# 4d create a dataset of only counties with annaul average PM2.5 >10 

# Summarizing Data
# 4e compute the correlation between aannual average PM2.5 and female unemployment
# 4f find the mean percent female unemployment among all sampled counties
# 4g compute the summary statistics of female unemployment
# 4h find the average percent female unemployment among the the sampled counties in NY 

#############################################
# hints: 
# in R, NA refers to missing data. 
# you can check whether a value is missing with is.na() 

a <- NA
is.na(a) 
a <- c(1,NA,2)
is.na(a)

# sometimes you need to tell functions to ignore or exclude missing values
# different functions have different arguments to address missing values 
# na.rm is common 
# often the ?function command can show you what arguments a function uses
# or you can check the full documentation

mean(df8$female.unemp)
mean(df8$female.unemp, na.rm= TRUE)
cor(df8$female.unemp, df8$male.unemp,use="complete")

# also, you can use objects inside of commands 

a   <- 100
df8 <- df8 %>% mutate(aveBMI_rescaled = aveBMI * a)

#############################################
#
#
###############################
#### 5: Plotting Data     #####
###############################

# Basic plotting functions.  

# 5a Histograms
# Histograms for county average BMI.  Describe the distribution of county mean BMI in this survey.  

?hist
hist(df8$aveBMI)

# we can set the number of bins with breaks option 

hist(df8$aveBMI, breaks=1)
hist(df8$aveBMI, breaks=20)

# we can set titles with xlab, ylab, and main 

hist(df8$aveBMI, breaks=20, 
     xlab = "Average BMI", ylab = "Count of Counties", main = "Histogram of Average BMI", las = 1) 

# las sets the axes labels orientation 
# we will cover how to arrange plots in the plotting workshop

# Instead of the histogram you can also plot density. 
# density is relative to the total data; the y axis is percentage not counts
# the density plots are smoothed over the range of data

plot(density(df8$aveBMI, na.rm = TRUE), col="red", las = 1)

# we can also overlay multiple plots on the same graph

hist(df8$aveBMI, breaks=20, xlab = "Average BMI", ylab = "Density", main = "Histogram of Average BMI", 
     freq = FALSE, las = 1) 
lines(density(df8$aveBMI, na.rm = TRUE), col = "red")

# 5b Boxplots
# Plot the distribution of annual PM among the sampled counties

boxplot(df8$avePM.idw, xlab = "Nationwide", ylab = expression("Annual PM"[2.5]))

# Plot the distribution of annual PM for each climatic region, using boxplot.  

boxplot(split(df8$avePM.idw, df8$climate.region), xlab = "NOAA Climate Region", ylab = expression("Annual PM"[2.5]))
boxplot(df8$avePM.idw ~ df8$climate.region,  xlab = "NOAA Climate Region", ylab = expression("Annual PM"[2.5]))

# 5c Scatterplots 
plot(df8$avePM.idw, df8$aveBMI, xlab = expression("Mean PM"[2.5]), ylab = "Average BMI", 
     main = "Scatterplot: \nCounty-average PM vs. BMI")

# we can add the least square fitted line to the above plot
# we will cover the lm() command in the next section

abline(lm(aveBMI ~ avePM.idw, data = df8), col = "blue", lwd = 2)

# 5d Multiple plots
# We can enhance the graph by specifying some options.  
# If you type help(plot) or help(par) for graphical parameters, you will find that there are many more 
# options that may be included that will allow you to adjust your graphs to your liking. 
# At times, you will find it convenient to have more than one graph on a single page.  
# R has a function that allows us to do that.  Here we divide the graphical page into one row and two 
# columns leaving two spaces for the following two plots.  
# It will remain setup this way until you reset it.

par(mfrow = c(1,2), las = 1)

boxplot(split(df8$avePM.idw, df8$climate.region), xlab = "NOAA Climate Region", ylab = expression("Annual PM"[2.5]),
        main = expression("PM"[2.5]*" Across Climate Regions"))

plot(df8$avePM.idw, df8$aveBMI, xlab = expression("Mean PM"[2.5]), ylab = "Average BMI", 
     main = "Scatterplot: County-average PM vs. BMI")

##################################
#### 6: Linear Regression    #####
##################################

# 6a Defining the question 

# Is county-average BMI associated with county-average annual PM2.5? 
# For this example I will use as outcome the average body mass index (BMI) which is a proxy for human body fat based on 
#     an individual's weight and height. Body mass index is defined as the individual's 
#     body weight divided by the square of his or her height therefore the units are kg/m2
# The county-average BMI was computed by taking a weighted average of survey responses in the
#     Behavioral Risk Factor Surveillance System (BRFSS), run by the CDC, where the weights are based on how 
#     representative the respondent is. 
# The PM data came from computing the inverse-distance-weighted average of nearby USEPA monitors

# 6b Initial look at data
# prior to creating the regression model, let's look at the distribution of aveBMI 

summary(df8$aveBMI)

# we can use the pdf command to save outputs of R into a pdf file

pdf(paste0(OutputPath, "Histograms.pdf"), width = 10) 

hist(df8$aveBMI, freq = FALSE, xlab = "BMI", ylab = "Density", las = 1)
lines(density(df8$aveBMI, na.rm = TRUE), col = "red")

dev.off()  # we use dev.off() to tell R to stop sending output to the pdf file. 

# we see that there is one county with very high aveBMI. 
# we can look at the data and see which observation has the highest BMI

LargestBMI <- df8 %>% filter(aveBMI == max(df8$aveBMI))

# 6c regression model 

mod <- lm(aveBMI ~ avePM.idw + per.black + per.latinx + per.asnam + 
            med.hinc + med.hval + lt.hs + female.unemp + male.unemp + climate.region, 
          data = df8, na.action = na.omit)

# notice that we used na.action to tell the function how to deal with na.omit 
# here lm() will ignore any observation that is missing data for any variable in the formula
# we will cover how to address missing data in Session 13: Missing Data

?lm

# 6d regression model results 
# we will first look at regression results and then consider modifying the model

summary(mod)

#   From the output you get: the model you run, estimates, SE and t-values of the parameters, 
#   and a description of the model in terms of residuals standard error, 
#   R-squared (proportion of total variation in the data which is explained by the model, 
#   and R-squared adjusted for the number of the parameter used in the fitted model) and F-test.

# what is in the model? what does the summary function produce? 

names(mod)
names(summary(mod))

# summary() and mod are lists, so we can extract specifc values

summary(mod)$coefficients
summary(mod)$coefficients[2,1] 
summary(mod)$coefficients[2,]  

# 6e Extracting results
# Sometimes it is useful to extract the whole object as a dataframe 

# we can write code to extract it ourselves
my.coeff <- summary(mod)$coefficients
# or
my.coeff <- data.frame(summary(mod)$coefficients)
names(my.coeff) <- c("est", "sterr", "tval", "pval")

# we can also use a packaged called broom to extract these values as a dataframe 
# install.packages("broom")
library(broom)

mod.tidy    <- tidy(mod)
modfit.tidy <- glance(mod)
confint_tidy(mod, conf.level = 0.95)
?confint_tidy
# note that broom does not work for every type of statistical model.

# 6f	Regression Diagnostics
#     Now let's see if the model follows the assumptions of a linear model. 
#     Four assumptions:
# a.	Linearity: constant slope
# b.	Normality
# c.	Independence
# d.	Constant variance
# 
# To check a): plot the data. We will learn how to fit more complex models such as using smoothing functions.
# To check b): look at the residuals. Use other plots such as the Q-Q plot (normal probability plots): 
#   the scatterplot should lie on the diagonal straight line.
# To check c) and d) plot of residuals vs fitted data. If the data are independent there should be no pattern 
#     in the data; if the variance is not constant you will see an increasing or decreasing cloud. 
#     We will learn what to do when these assumptions are violated.
#
# For this model, we will also remove outliers that have very large residuals or very high leverage

# In the case of linear model, the plot of the model gives diagnostic plots

par(mfrow = c(2,2))
plot(mod)
# The same plots can be obtained by the following:

par(mfrow=c(3,1))
plot(mod$residuals)
hist(mod$residuals)
plot(density(mod$residuals))

## get the residuals and plot against fitted
par(mfrow=c(2,2))
res.mod    <- residuals(mod) # or res.bp.model    <- bp.model5$residuals
fitted.mod <- fitted(mod)    # or fitted.bp.model <- bp.model5$fitted.values
plot(fitted.mod, res.mod)

## QQ plot
qqnorm(res.mod)
qqline(res.mod)

## get the square root standardized residuals and plot against fitted
res.sqsr <- sqrt(abs(scale(resid(mod))))
plot(fitted.mod, res.sqsr)

## get the leverage and plot against standardized residuals 
lev    <- hat(model.matrix(mod))
res.sr <- scale(resid(mod))
plot(lev, res.sr, las = 1, ylab = "Standardized Residuals", xlab = "Leverage")
abline(h = 0, col = "red")

# Let's see if our model meets the assumptions
par(mfrow=c(2,2))
plot(mod)
# a: We will cover this in session 3: NonLinearity
# b: We can look at this guide to help us interpret the QQ plot 
#    http://seankross.com/2016/02/29/A-Q-Q-Plot-Dissection-Kit.html
#    It appears that the tails of the distribution of our residuals is a bit fatter 
#    than would be expected from a normal distribution
# c: The residuals do not appear to increase or decrease in the residuals vs fitted plot
# d: The variance of residuals appears consistent across predicted values of aveBMI
# outliers: 
# based on the plot() results, we see that observation #154 have a leverage of 1, so we will see 
# what happens when we remove it 
# from the residuals vs fitted, we see that observation #169 has a very large residual, we'll also
# check how influential this point is

# 6f Modified regression 

# look at observation 230 
df8[169,]

# remove the potential outliers
df9 <- df8[-c(169, 154),]

# the regression model 
mod1 <- lm(aveBMI ~ avePM.idw + per.black + per.latinx + per.asnam +  
             med.hinc + med.hval + lt.hs + female.unemp + male.unemp + climate.region, 
           data = df9, na.action = na.omit)


# compare models 
# did removing the outliers change the model? 

summary(mod1)
summary(mod)

# Yes, several coefficients changed by a fair amount 
# so we will use df9, without the outliers. 
# removing outliers should not be done automatically, we need to think about it. 
# and we should never decide to remove outliers just because it changes the model to agree with our hypothesis!

# 6g Interaction term 
# Does the association between county-level annual PM2.5 and county average BMI vary by 
# climate region? 
# Let's add an interaction term 

# regression model 

mod2 <- lm(aveBMI ~ avePM.idw*climate.region + per.black + per.latinx + per.asnam +  
             med.hinc + med.hval + lt.hs + female.unemp + male.unemp , 
           data = df9, na.action = na.omit)

# Let's look at the model 

summary(mod2)

# we see that the term for pm*ohio_valley is almost significant 
# We can use the values from the model to estimate 
# the PM-aveBMI association within the Ohio Valley
# and to estimate the confidence intervals of the association 

# 6h Compute the association and its uncertainty 

# here we create tables of coefficients and covariance
coef.mat <- summary(mod2)$coefficients
var.mat  <- vcov(mod2)

# the total term for the association is the 
# sum of the term in the reference region plus the term for Ohio Valley

beta.ohio_valley <- coef.mat["avePM.idw",1] + 
  coef.mat["avePM.idw:climate.regionohio_valley",1]

# Compute variance in order to compute standard error
# We must compute the variance for the total term 
# Var(Beta1 + Beta3) = Var(Beta1) + Var(Beta3) + CoVar(Beta1, Beta3) + CoVar(Beta3, Beta1)
# Var(Beta1 + Beta3) = Var(Beta1) + Var(Beta3) + 2*CoVar(Beta1, Beta3) 

var.ohio_valley  <- var.mat["avePM.idw", "avePM.idw"] + 
  var.mat["avePM.idw:climate.regionohio_valley", "avePM.idw:climate.regionohio_valley"] +
  2*var.mat["avePM.idw", "avePM.idw:climate.regionohio_valley"]

ste.ohio_valley  <- sqrt(abs(var.ohio_valley))

# compute confidence intervals 

lci.ohio_valley <- beta.ohio_valley - 1.96*ste.ohio_valley
uci.ohio_valley <- beta.ohio_valley + 1.96*ste.ohio_valley

OHval <- paste(round(beta.ohio_valley, 3), " (95% CI: ", round(lci.ohio_valley, 3), ", ", 
               round(uci.ohio_valley, 3), ")", sep = "")

# we can ask whether including the interaction terms improved model fit 
# anova() provides a nice test 

anova(mod1, mod2)

# the results of the anova are not statistically significant, 
# indicating that as a whole interaction between pm and region does not improve model fit
# i.e. no evidence of effect modification by region

# 6i Stratification 
# we can stratify with the subset command 

mod3 <- lm(aveBMI ~ avePM.idw + per.black + per.latinx + per.asnam +  
             med.hinc + med.hval + lt.hs + female.unemp + male.unemp , 
           data = df9, na.action = na.omit, subset = (climate.region == "ohio_valley"))

summary(mod3)

###############################
#### 7: Exercises         #####
###############################

# Plotting 
# 7a Create a histogram of median household value
# see if you can make the bars orange, or a color of your choice: 
# http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

# 7b Create a boxplot of average BMI within levels of median home value 
# hint: first make a new variable of 4 quantiles of median home value 

# 7c Create a scatterplot of median home value and average bmi 
# add a trend line 
# see if you can make the trend line orange

# Modeling 

# 7d.i Compute the confidence intervals for avePM.idw in the model 
#       aveBMI~ avePM.idw + med.hval +female.unemp + male.unemp
# 7d.ii Compute the term and confidence intervals for a 10 unit increase in avePM.idw

# 7e According to mod2, compute the term and confidence intervals for avePm.idw in the southeast region

# 7f Does the coefficient for avePM.idw change when we do not control for the potential confounders? 

# 7g Is the relationship between annual average PM and average BMI the same in 
#    counties with below average male unemployment, and those with above average male unemployment? 

#############################################
# hints: 
# when we multiply a term by a unit, we also multiply the standard errors 
# for example, we use 
# 10*(beta.ohio_valley + 1.96*ste.ohio_valley) 
# not ( 10*beta.ohio_valley + 1.96*ste.ohio_valley)

# what does an unadjusted model look like? 

# you will need to make a categorical variable for male unemployment. 

#############################################
####################################################################################

# optional questions (if time) 

# 7e.ii Compute the term and confidence intervals for a 10 unit increase in avePM.idw in the southeast region

# 7e Test whether adding an interaction term for education (lt.hs) and PM2.5 improves the model 

# 7f Assuming the interaction between lt.hs and PM2.5, compute the main effect and confidence intervals
# for a 10 unit increase of pm, for a county with 15% less than Hs education 

