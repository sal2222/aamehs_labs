# CaseCrossover Session 
# 1: Cleaning the data 
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# Feburary 8, 2019
# Session 4: CaseCrossover
##################################################
#### Table of Contents ####
# 0: Preparation 
# 1: Prepare Data 
# 2: Combine into Daily Observations
# 3: Organize Data
# 4: Create Control Days 
# 5: Assign Exposure to Days
# 6: Finalize Data
########################
#### 0: Preparation ####
########################
# 0a Install packages 

# mgcv allows for quick modeling of generalized additive models 
# install.packages("mgcv")

# 0b Load packages
library(readr)
library(dplyr) 
library(ggplot2)
library(lubridate)
library(stringr)
library(tidyr)

#0c Declare directories

DataPath   <- "~/Desktop/AAMEHS/session4_CaseCrossover/Lab/Data/"
OutputPath <- "~/Desktop/AAMEHS/session4_CaseCrossover/Lab/Outputs/"

# 
# DataPath   <- "/Users/marianthi_anna/Dropbox/AAMEHS_Spring2019/Labs/Session2_QuantileRegression/"
# 
# 
###############################################################
##########################
##### 1: Prepare Data ####
##########################
# 1a Readin Data
df <- read_csv(paste0(DataPath, "processed_nyc_aqs_cvd_mh.csv"), 
               col_types = cols(.default ="c"))

# 1b Keep only the variables we want 
df <- df %>% select(county, dateD, pm_aqs, temp, rh, cvd_in) 

# 1c Keep only NYC 
NYCcounties <- c("Bronx", "Kings", "New York", "Queens", "Richmond")
df <- df %>% filter(county %in% NYCcounties)

# 1d Convert to numeric 
df <- df %>% mutate(pm_aqs = as.numeric(pm_aqs),
                      temp = as.numeric(temp), 
                      rh = as.numeric(rh),
                      cvd_in = as.numeric(cvd_in))

#############################################
##### 2: Combine into Daily Observations ####
#############################################
# 2a add 2010 burough population
# data from https://www1.nyc.gov/site/planning/data-maps/nyc-population/current-future-populations.page

# create burough dataframe
totpop <- 8175133
burough <- c( "Bronx", "Brooklyn", "Manhattan", "Queens", "Staten Island")
county <- c("Bronx", "Kings", "New York", "Queens", "Richmond")
pop <- c(1385108, 2504700, 1585873, 2230722, 468730)

countypop <- data.frame(burough, county, pop, stringsAsFactors = FALSE)

# 2b Compute population weights

countypop$weight <- countypop$pop/totpop

# 2c Assign population weights 

df1 <- df %>% left_join(countypop, by = "county")

# 2d Combine observations from buroughs, weighted by population 

nyc.df <- df1 %>% group_by(dateD) %>% 
  summarize(dailyPM = sum(pm_aqs*weight), 
            dailyTemp = sum(temp*weight),
            dailyRH = sum(rh*weight), 
            dailyCVDin = sum(cvd_in))

##########################
#### 3: Organize Data ####
##########################

# 3a Remove days with missing data
nyc.df <- nyc.df %>% filter(complete.cases(nyc.df))

# 3b Separate the population and exposure data
pop <- nyc.df %>% select(dailyCVDin, dateD)
exp <- nyc.df %>% select(dateD, dailyPM, dailyTemp, dailyRH)

# 3c Extract the date values for the exp 
exp1 <- exp %>% mutate(ExpDate = parse_date_time(dateD, "mdy")) %>% select(-dateD)
pop1 <- pop %>% mutate(CaseDate = parse_date_time(dateD, "mdy")) %>% select(-dateD)

################################
#### 4: Create Control Days ####
################################

# We need to identify the control days for each case day,
# and then we will assign exposures for each day 

# Creating Control Days 
# we first create a set of potential control days 
# that are up to 5 weeks before and after the case day 
# the potential control days are the same day of the week 
# as the case day 

# 4a Define function to create control days 
make_control_day <- function(BeforeAfter, WK){ 
  VarName <- paste0(BeforeAfter, "_", str_trunc(WK,1,"left",""))          #the name of the lag 
  pop1 %>% mutate(!!VarName := CaseDate + as.period(7*WK, "day"))   #adds WKs number of weeks, preserves hour of day even in daylightsaving
}

# 4b Create days 
pop1 <- make_control_day("Before", -4)
pop1 <- make_control_day("Before", -3)
pop1 <- make_control_day("Before", -2)
pop1 <- make_control_day("Before", -1)
pop1 <- make_control_day("CaseDate", 0)
pop1 <- make_control_day("After", 1)
pop1 <- make_control_day("After", 2)
pop1 <- make_control_day("After", 3)
pop1 <- make_control_day("After", 4)

# 4c Put in long format by DayName
days <- pop1 %>% gather("DayName", "DayDate", contains("CaseDate_"),contains("Before_"), contains("After_") ) 

# 4d Stratify by month of event 
days1 <- days %>% filter(month(CaseDate) == month(DayDate))

#assign exposure to each day... 

#####################################
#### 5: Assign Exposure to Days  ####
#####################################

# 5a Assign exposure based on date
days_exp <- days1 %>%
  inner_join(exp1, by = c("DayDate" = "ExpDate"))

##########################
#### 6: Finalize Data ####
##########################
# 6a Add variables for case-crossover analysis 
# case =1 if it is a case
days_exp <- days_exp %>% 
  mutate(  Case = if_else(str_detect(DayName, "Case"), 1, 0),
           TimetoEvent = if_else(str_detect(DayName, "Case"), 1, 2), 
           DayID = as.character(CaseDate))

# 6b Select the variables we want to keep 

NYC_daily <- days_exp %>% select(DayID, DayName, DayDate, Case, TimetoEvent, dailyPM, dailyTemp, dailyRH, dailyCVDin)

# 6c Remove missing years
# there is a lot of missingness in 2011 and 2012, so we will drop those years
# evidence: 
ggplot(NYC_daily, aes(DayDate, dailyTemp)) + 
  geom_line()

# drop years 
NYC_daily <- NYC_daily %>% filter(year(DayDate) < 2010)

# 6d Save results

NYC_daily %>% write_csv(paste0(DataPath, "nyc_daily_processed.csv"))


