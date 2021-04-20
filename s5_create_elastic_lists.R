## Elastic Net Regression for MOTIVC resting-state fMRI data
#Using formula for: Varying signals. High correlation between predictors
setwd("E:\\0_parcellation_analysis")
library("ff")
library(pacman)
p_load(MASS)  # Package needed to generate correlated precictors
p_load(glmnet)  # Package to fit ridge/lasso/elastic net models


#load csv files with difference data between dyads for every node-to-node correlation
#the size of these files are as follows:
# lower 4 2017 girls = 254x36857 double
# upper 4 2017 girls = 379x36857 double
# lower 4 2018 girls = 137x36857 double 
# these sizes apply to both distance and community similarity files
##these were last updated using mutual SN information on June 5 2019

years <- c("Cohort 1", "Cohort 3", "Cohort 2")
filenames <- c("MODELLING_L4.csv","MODELLING_4.csv","MODELLING_700.csv")
data_nolabels <- lapply(X = filenames, FUN = read.csv, header = F)
names(data_nolabels) <- years

filenames <- c("yL4_dist.csv","y4_dist.csv","y700_dist.csv")
dist_labels <- lapply(X = filenames, FUN = read.csv, header = F)
names(dist_labels) <- years

filenames <- c("yL4_comm.csv","y4_comm.csv","y700_comm.csv")
comm_labels <- lapply(X = filenames, FUN = read.csv, header = F)
names(comm_labels) <- years

save(list=c("data_nolabels","dist_labels", "comm_labels"),file = 'shiny_for_elastic.RData')


#ONCE YOU HAVE CREATED THE RData FILE ABOVE, YOU CAN JUST COMMENT OUT THE ABOVE SCRIPT 
#AND LOAD THE FILE INSTEAD BY UNCOMMENTING THE NEXT LINE
#load("shiny_for_elastic.RData")

#CREATE LIST FOR COMMUNITY DATA
community_data <- list()


for (yr in years){
# Generate data
set.seed(19873)
n <- nrow(data_nolabels[[yr]])    # Number of observations
p <- ncol(data_nolabels[[yr]])     # Number of predictors included in model

# Split data into train and test sets
# alternatively, you can use one year for training and another year for testing
train_rows <- sample(1:n, .66*n)
community_data[[yr]]$x.train <- as.matrix(data_nolabels[[yr]][train_rows, ])
community_data[[yr]]$x.test <- as.matrix(data_nolabels[[yr]][-train_rows, ])

#decide on whether you want community similarity or social distance metrics used for analysis
community_data[[yr]]$y.train <- as.factor(comm_labels[[yr]][train_rows,1])
community_data[[yr]]$y.test <- as.factor(comm_labels[[yr]][-train_rows,1])
# y.train <- dist_labels[[dat]][train_rows,1]
# y.test <- dist_labels[[dat]][-train_rows,1]
}
distance_data <- list()
for (yr in years){
  # Generate data
  set.seed(19873)
  n <- nrow(data_nolabels[[yr]])    # Number of observations
  p <- ncol(data_nolabels[[yr]])     # Number of predictors included in model
  
  # Split data into train and test sets
  train_rows <- sample(1:n, .66*n)
  distance_data[[yr]]$x.train <- as.matrix(data_nolabels[[yr]][train_rows, ])
  distance_data[[yr]]$x.test <- as.matrix(data_nolabels[[yr]][-train_rows, ])
  
  #decide on whether you want community similarity or social distance metrics used for analysis
  # y.train <- as.factor(comm_labels[[dat]][train_rows,1])
  # y.test <- as.factor(comm_labels[[dat]][-train_rows,1])
  distance_data[[yr]]$y.train <- dist_labels[[yr]][train_rows,1]
  distance_data[[yr]]$y.test <- dist_labels[[yr]][-train_rows,1]
}
## Fitting models

fits_community <- list()
for (yr in years){
## if using community structure, use binomial settings
fits_community[[yr]]$fit.lasso <- glmnet(community_data[[yr]]$x.train, community_data[[yr]]$y.train, family="binomial", alpha=1)
fits_community[[yr]]$fit.ridge <- glmnet(community_data[[yr]]$x.train, community_data[[yr]]$y.train, family="binomial", alpha=0)
fits_community[[yr]]$fit.elnet <- glmnet(community_data[[yr]]$x.train, community_data[[yr]]$y.train, family="binomial", alpha=.5)
# 10-fold Cross validation for each alpha = 0, 0.1, ... , 0.9, 1.0
fits_community[[yr]]$fit.lasso.cv <- cv.glmnet(community_data[[yr]]$x.train, community_data[[yr]]$y.train, type.measure="mse", alpha=1,
                           family="binomial")
fits_community[[yr]]$fit.ridge.cv <- cv.glmnet(community_data[[yr]]$x.train, community_data[[yr]]$y.train, type.measure="mse", alpha=0,
                          family="binomial")
fits_community[[yr]]$fit.elnet.cv <- cv.glmnet(community_data[[yr]]$x.train, community_data[[yr]]$y.train, type.measure="mse", alpha=.5,
                          family="binomial")

for (i in 0:10) {
    assign(paste("comm_fit",yr, i, sep="_"), cv.glmnet(community_data[[yr]]$x.train, community_data[[yr]]$y.train, type.measure="mse",
                                              alpha=i/10,family="binomial"))}
}

fits_distance <- list()
for (yr in years){
## if using social distance, use gaussian settings
fits_distance[[yr]]$fit.lasso <- glmnet(distance_data[[yr]]$x.train, distance_data[[yr]]$y.train, family="gaussian", alpha=1)
fits_distance[[yr]]$fit.ridge <- glmnet(distance_data[[yr]]$x.train, distance_data[[yr]]$y.train, family="gaussian", alpha=0)
fits_distance[[yr]]$fit.elnet <- glmnet(distance_data[[yr]]$x.train, distance_data[[yr]]$y.train, family="gaussian", alpha=.5)
fits_distance[[yr]]$fit.lasso.cv <- cv.glmnet(distance_data[[yr]]$x.train, distance_data[[yr]]$y.train, type.measure="mse", alpha=1,
                          family="gaussian")
fits_distance[[yr]]$fit.ridge.cv <- cv.glmnet(distance_data[[yr]]$x.train, distance_data[[yr]]$y.train, type.measure="mse", alpha=0,
                          family="gaussian")
fits_distance[[yr]]$fit.elnet.cv <- cv.glmnet(distance_data[[yr]]$x.train, distance_data[[yr]]$y.train, type.measure="mse", alpha=.5,
                          family="gaussian")

for (i in 0:10) {
  assign(paste("dist_fit",yr, i, sep="_"), cv.glmnet(distance_data[[yr]]$x.train, distance_data[[yr]]$y.train, type.measure="mse",
                                            alpha=i/10,family="gaussian"))
}
}

fits <- list(fits_community,fits_distance)
data <- list(community_data,distance_data)
names(fits) <- c("community","distance")
names(data) <- c("community","distance")
save(list=c("fits", "data"),file = 'ElasticNetReg.RData')

#ONCE YOU HAVE CREATED THE ABOVE RDATA FILE, YOU CAN COMMENT OUT ALL THE ABOVE AND JUST LOAD 
#THE DATA INSTEAD BY UNCOMMENTING THE NEXT LINE
#load('ElasticNetReg.RData')


## Plot solution path and cross-validated MSE as function of lambda

type <- c("distance")
yr <- c("Cohort 3")
# Plot solution paths:

par(mfrow=c(3,2))
# For plotting options, type '?plot.glmnet' in R console
plot(fits[[type]][[yr]]$fit.lasso, xvar="lambda")
plot(fits[[type]][[yr]]$fit10, main="LASSO")

plot(fits[[type]][[yr]]$fit.ridge, xvar="lambda")
plot(fits[[type]][[yr]]$fit0, main="Ridge")

plot(fits[[type]][[yr]]$fit.elnet, xvar="lambda")
plot(fits[[type]][[yr]]$fit5, main="Elastic Net")





## MSE on test set


yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test)
yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test)
yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test)
yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test)
yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)
yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test)
yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test)
yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test)
yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)


#the following tests will only be relevant for gaussian models
mse0 <- mean((y.test - yhat0)^2)
mse1 <- mean((y.test - yhat1)^2)
mse2 <- mean((y.test - yhat2)^2)
mse3 <- mean((y.test - yhat3)^2)
mse4 <- mean((y.test - yhat4)^2)
mse5 <- mean((y.test - yhat5)^2)
mse6 <- mean((y.test - yhat6)^2)
mse7 <- mean((y.test - yhat7)^2)
mse8 <- mean((y.test - yhat8)^2)
mse9 <- mean((y.test - yhat9)^2)
mse10 <- mean((y.test - yhat10)^2)



