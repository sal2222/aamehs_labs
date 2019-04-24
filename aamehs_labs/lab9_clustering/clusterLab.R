# K-Means and Hierarchial Clustering Lab Session
# Marianthi-Anna Kioumourtzoglou
# Sebastian Rowland
# Advanced Analytic Methods for Environmental Epidemiology
# March 27, 2019
# Session 9: K-means Hierarchial Clustering 
# Special thanks to 
# Yanelli Nunez, Lizzy Gibson, Ahlam Abuawad, and Marianthi Kioumurtzoglou 
# For example code from the 2018 Mixtures Workshop

####***********************
#### Table of Contents ####
# 0:  Preparation 
# 1:  Prepare Exposure Data

# A: K-Means
# 2: K-Means Solution
# 3: Sum of Squares 
# 4: Plot K-Means Solution 
# 5: Create Exposure Matrix of K-Means Clusters

# B: Hierarchial Clustering
# 6: Hierarchial Cluster Solution 
# 7: Choose Optimal Tree Depth 
# 8: Plot Hierarchial Cluster Solution
# 9: Create Exposure Matrix of Hierarchial Clusters

####********************
#### 0: Preparation ####
####********************
# 0a Install packages 

# install.packages("reshape2")
install.packages("dendextend")
install.packages("ggdendro")

# 0b Load packages
library(lubridate)
library(tidyverse)

library(dendextend)
library(ggdendro)
library(factoextra)
library(reshape2)

####*******************************
#### 1: Prepare Exposure Data #####
####*******************************

# 1a Readin Data 

df0 <- read_csv("Boston_PM_constituents.csv")

# 1b Convert date variable to datetime format 

df0$Date <- parse_date_time(df0$Date, "mdy")

# 1c View data 

head(df0)

# 1e.i Create a list of days to remove

FireWorksDays <- c("7_2", "7_3", "7_4", "7_5")

# 1e.ii Remove those days, for each year 

df0 <- df0 %>% mutate(MonthDay = paste0(month(Date), "_", day(Date))) %>% 
  filter(!(MonthDay %in% FireWorksDays )) %>% 
  dplyr::select(-MonthDay)

# 1f Remove Variables
# we will also remove Sodium
# because XRF readings for Na are not so accurate 

df0 <- df0 %>% dplyr::select(-Sodium)

# 1g Remove days with missing data 

df0 <- df0 %>% filter(complete.cases(df0))

df <- df0 %>% dplyr::select(-Date, -ent)

# 1i Check dimensions of dataset

dim(df)

# 1j Summary statistics on dataset

summary(df)

# 1k Scale the dataset 

df <- as.data.frame(scale(df, center = TRUE))

# 1l Correlation Plot

cormat <- cor(df, use = "pairwise.complete.obs", method = c("spearman"))

melted_cormat <- melt(cormat) %>% dplyr::rename(Correlation = value)

melted_cormat <- melted_cormat  %>% 
  ggplot(aes(Var1, Var2, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-0.1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation\n") +
  theme_minimal(base_size = 15) + # minimal theme
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, size = 12, hjust = 1),
        axis.text.y = element_text(vjust = 0.3, size = 12, hjust = 1),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15)) +
  coord_fixed() + labs(x = "", y = "")


melted_cormat

####****************
#### A: K-Means ####
####****************


####*************************
#### 2: K-means Solution ####
####*************************


# 2a Set seed  
# k-means places the initial centorids at random locations, 
# if we set the seed, then each of our versions of R will 
# initialize with the same sets of starting points 

set.seed(42)

# 2b K-means with 2 clusters 
# the solution of k-means can depend on the initial locations
# of the starting points 
# so we will try multiple initial locations 
# and choose the best solution 
# The best solution is defined as the solution with the least 
# total within-cluster sum of square error 
# (between the cluster mean and each member of the cluster)


km.2 <- kmeans(df,
                centers = 2,   # number of clusters 
                nstart  = 100, # number of starts for random assignments
                iter.max = 20) # number of iterations

# 2c Features of K-means Solution  

# size of each cluster 
km.2$size

# cluster assignment for each observation
head(km.2$cluster)

# Mean concentration of each constituent 
# within each cluster
# note that these have been scaled!
km.2$centers

# total sum of squares
km.2$totss 

# within-cluster sum of squares
# for each cluster
km.2$withinss 

# Total within-cluster sum of squares
# this is what we want to minimize! 
km.2$tot.withinss 
sum(km.2$withinss)

# Between-cluster sum of squares
km.2$betweenss 

# Percentage of total variance that comes from between-cluster variance
Percent_bt_k2 = paste(round(100*(km.2$betweenss/km.2$totss),1), "%", sep = "")
Percent_bt_k2

# 2d K-means with 6 clusters 

km.6 <- kmeans(df,
                centers = 6,   # number of clusters 
                nstart  = 100, # number of random locations to assess
                iter.max = 20) 

# 2e Features of K-means solution
km.6$size
head(km.6$cluster) ## cluster assignment
km.6$centers ## center means
km.6$totss ## total sum of squares
km.6$withinss ## within cluster sum of squares by cluster
km.6$tot.withinss ## total within cluster sum of squares -- this is what we want to minimize!
km.6$betweenss ## between cluster sum of squares
Percent_bt_k6 = paste(round(100*(km.6$betweenss/km.6$totss),1), "%", sep = "")


# 2f Compute variance for a range of k-means solutions

# use a loop to run k-means for 1 - 50 clusters and generate variance values

km.res <- data.frame(matrix(NA, 50, 3))

# warning: this takes a moment to run 

for (k in 1:50) {
  kk <- kmeans(df, k, nstart = 50, iter.max = 20)
  km.res[k, ] <- cbind(k, kk$tot.withinss, kk$totss) 
}

# Name the variables created
names(km.res) <- c("clusters", "WithinSS", "TotSS")

# Create a variable for within cluster variance/total SS
km.res$PropWithin <- 100*km.res$WithinSS/km.res$TotSS

# Look at table of variances
km.res

# Plot a subset of data
km.res  %>% 
  ggplot(aes(x = clusters, y = PropWithin)) + geom_line() + 
  geom_vline(xintercept = 6, color = "red", linetype = "dotted") + 
  labs(title = "", 
       y = "Proportion of Within over Total SS", 
       x = "Number of Clusters")


####################################
####*******************************
#### 4: Plot K-Means Solution  ####
####*******************************

# 4a Create dataframe of means of clusters 

# 4a.i Extract cluster means from solution

km6_centers <- as.data.frame(t(km.6$centers))

# 4a.ii Create column of species names

km6_centers$species <- row.names(km6_centers)

# 4a.iii Put data in long format

plot_means_km6 <- km6_centers %>% 
                   gather(key = "Cluster", value = "mean", -species) 

# 4a.iv Rename clusters

plot_means_km6 <- plot_means_km6 %>% 
                   mutate(Cluster = fct_recode(Cluster, 
                              "Cluster 1" = "1", 
                              "Cluster 2" = "2", 
                              "Cluster 3" = "3", 
                              "Cluster 4" = "4", 
                              "Cluster 5" = "5", 
                              "Cluster 6" = "6"))

# 4.b Plot means of each cluster

ggplot(plot_means_km6, aes(x = species, y = mean, fill = species)) +
  geom_col() +
  geom_hline(yintercept = 0, size = 0.2) +
  facet_wrap(~ Cluster) +        # creates 6 plots, for each cluster
  theme_bw() 

##*****************************************************
#### 5: Create Exposure Matrix of K-Means Clusters #####
####****************************************************

# We can combine the date variable and the cluster assignment to 
# create a new exposure matrix based on which cluster the day belongs to 
# We can then use this exposure matrix in an epidemiological analysis 
# where each cluster is treated as a categorical dummy variable. 

# 5a Extract cluster membership 

clusters.df <- as.data.frame(km.6$cluster)

colnames(clusters.df) <- "cluster"

# 5b Convert clusters to character 
# unlike in PCA, the number of the clusters do not indicate order 
# they are just unique identifiers 

clusters.df$cluster <- as.factor(clusters.df$cluster)

# 5c Add the dates from the original data 

clusters.df$Date <- df0$Date

# 5d View dataframe 

head(clusters.df)

# 5e Look for weather patterns in clusters 

# 5e.i Read in weather data 

weather <- read_csv("noaa_boston_weather.csv")

# 5e.ii Organize the date column 

weather <- weather %>% dplyr::rename(Date = DATE)
weather <- weather %>% mutate(Date = parse_date_time(Date, "ymd"))

# 5e.iii Join weather data and cluster data 

df.wea <- clusters.df %>% 
          left_join(weather, by = "Date")

# 5e.iv Compute mean weather for each cluster 

summary(df.wea)

df.wea.means <- df.wea %>% 
                group_by(cluster) %>% 
                summarize(
                meanWind   = mean(AWND), 
                meanPrecip = mean(PRCP, na.rm = TRUE), 
                meanSnow   = mean(SNOW, na.rm = TRUE), 
                meanTMin   = mean(TMIN), 
                meanTMax   = mean(TMAX))

# 5e.v  Put data in long format

df.wea.means.long <- df.wea.means %>% 
                     gather(WeatherVar, "Value", -cluster)

# 5e.vi Plot average weather of each cluster 

df.wea.means.long %>% 
  ggplot(aes( cluster, Value)) + 
  geom_point(aes(color = WeatherVar), shape = 18, size = 3) + 
  facet_wrap(~WeatherVar, scales = "free")

####*******************************
#### B: Hierarchial Clustering ####
####*******************************

# clean the environment 

rm(km.6, km6_centers, plot_means_km6, Percent_bt_k6, clusters.df)


####*************************************
#### 6: Hierarchial Cluster Solution ####
####*************************************


df.dist <- dist(df, method = "euclidean")

# 6a.ii Create hierarchial cluster solution
# hclust is an agglomerative clustering algorithm
# other functions are available in R 
# the method arguement determines the function used to compute distance between clusters

hc.complete <- hclust(df.dist, method = "complete")

# 6b Look at dendrogram

as.dendrogram(hc.complete) %>% head()

####******************************************
#### 7: Plot Hierarchial Cluster Solution ####
####******************************************
# hierarchial clusters are most clearly represented as visual trees

# 7a Plotting with Base R 

# 7a.i Extract the dendrogram from the HC solution 

dendro.complete <- as.dendrogram(hc.complete)

# 7a.ii Plot 
# height indicates how dissimilar the two clusters are 
# i.e. the clusters fused at height = 22 
# are much more dissimilar 
# than the clusters fused at height = 4

dendro.complete %>% 
  plot(main = "Complete Linkage", ylab = "Height", leaflab = "none") 


# 7b Plotting with GGplot

# 7b.i Extract branch node and height data as a dataframe
# by extracting these as a dataframe, and not a dendrogram object 
# we can input it into ggplot 

# this line will take a minute to run

# dendro.complete.df <- dendro_data(hc.complete, type = "rectangle")

# 7b.ii Plot

# plot.complete <- ggplot(segment(dendro.complete.df)) + 
#  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), alpha = 0.5) +
#  theme_minimal() +
#  theme(axis.text.x = element_blank(), 
#        axis.ticks.x = element_blank(), 
#        panel.grid = element_blank()) + 
#  labs(y = "Height", x = "", title = "Complete Linkage") 
#
# plot.complete

####***********************************
#### 8: Choose Optimal Tree Depth  ####
####***********************************

# 8a Plot branches of different cut heights
# we can use the color_branches() function to divide our data 
# into various numbers of clusters
# the dashed line (which we manually create with abline())
# indicates where the cut is taking place

dendro.complete  %>% 
  color_branches(k = 2) %>%
  plot(main = "Complete Linkage", ylab = "Height", leaflab = "none") 
abline(h = 22, lty = 3)

dendro.complete  %>% 
  color_branches(k = 7) %>%
  plot(main = "Complete Linkage", ylab = "Height", leaflab = "none") 
abline(h = 16, lty = 3)

dendro.complete  %>% 
  color_branches(k = 20) %>%
  plot(main = "Complete Linkage", ylab = "Height", leaflab = "none") 
abline(h = 9.3, lty = 3)

# 8b Plot rectangles around clusters 

# 8b.i Plot with 4 clusters

dendro.complete %>% 
  plot(main = "Complete Linkage", ylab = "Height", leaflab = "none") 

rect.hclust(hc.complete, k = 4, border = 2:5)

# 8b.ii Plot with 12 clusters

dendro.complete %>% 
  plot(main = "Complete Linkage", ylab = "Height", leaflab = "none") 

rect.hclust(hc.complete, k = 12, border = 2:5)

# 8c Plot mean values for each cluster

# 8c.i Cut the tree to have 12 clusters

hc.cluster.12 <- cutree(hc.complete, 12)
table(hc.cluster.12)

# 8c.ii Include the cluster assignments in the original scaled data

df.hc <- df %>% mutate(hcluster.12 = hc.cluster.12)

# 8c.iv Compute mean concentration of each constituent within each cluster

df.hc.mean <- df.hc %>% 
  group_by(hcluster.12) %>% 
  summarize_all(.funs = mean)

# 8c.v Put data in long format

plot_means_hc12 <- df.hc.mean %>%
  gather(key = "species", value = "mean", -hcluster.12) 

# 8c.vii Plot means of each cluster
# we will just plot 7 of the clusters for aesthetic reasons 

# 8c.vii.i Randomly select clusters to view

plot_means_hc12.7c <- plot_means_hc12 %>% 
  filter(hcluster.12 %in% 1:6)

# 8c.vii.ii Plot

ggplot(plot_means_hc12.7c, aes(x = species, y = mean, fill = species)) +
  geom_col() +
  geom_hline(yintercept = 0, size = 0.2) +
  facet_wrap(~ hcluster.12) +        # creates 12 plots, for each cluster
  theme_bw() 

###*********************************************************
#### 9: Create Exposure Matrix of Hierarchial Clusters #####
####********************************************************

# We can combine the date variable and the cluster assignment to 
# create a new exposure matrix based on which cluster the day belongs to 
# We can then use this exposure matrix in an epidemiological analysis 
# where each cluster is treated as a categorical dummy variable. 

# 9a Extract cluster membership 

hc.cluster.12.df <- as.data.frame(cutree(hc.complete, 12))

colnames(hc.cluster.12.df) <- "cluster"

# 9b Convert clusters to character 
# unlike in PCA, the number of the clusters do not indicate order 
# they are just unique identifiers 
# in fact, the numbers will vary each time you run kmeans()
# We will convert the clusters to letters to make sure that 
# R will treat them as dummy variables 
# and so that anyone else who looks at our data knows that they are categories 
# if you have more that 26 clusters, you will need to write a different function 

# 9c Add the dates from the original data 

hc.cluster.12.df$Date <- df0$Date

# 9d View dataframe 

head(hc.cluster.12.df)



