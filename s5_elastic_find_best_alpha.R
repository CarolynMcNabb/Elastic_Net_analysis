## Carolyn McNabb
##Elastic Net Regression for MOTIVC resting-state fMRI data
#Using formula for: Varying signals. High correlation between predictors
home_dir <- c("F:\\0_parcellation_analysis/scripts-data-sharing/")
setwd(home_dir)
library(pacman)
p_load(MASS)  # Package needed to generate correlated precictors
p_load(glmnet)  # Package to fit ridge/lasso/elastic net models
p_load("ff")
#if you need to install metaforlmer, uncomment the two lines below
# library(devtools)
# install_github("https://github.com/LilyFG/metaforlmer.git")
library(lme4)
library(metaforlmer)
library(metafor)
library(lmerTest)
#library(ggforest)
# library(survminer)
#p_load(SDMTools)
# require(doMC)
# registerDoMC(cores=2)
#load csv files with difference data between dyads for every node-to-node correlation
#the size of these files are as follows:
# lower 4 2017 girls = 254x36857 double
# upper 4 2017 girls = 379x36857 double
# lower 4 2018 girls = 137x36857 double 
# these sizes apply to both distance and community similarity files
##these were last updated using mutual SN information on June 5 2019
years <- c("Cohort 1", "Cohort 3", "Cohort 2")

#YOU ONLY NEED TO RUN THE BELOW CODE ONCE AND THEN YOU CAN COMMENT OUT AND 
#INSTEAD LOAD THE SAVED DATA

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
load("shiny_lists.RData")

#########################################
#model a
x.train <- as.matrix(rbind(data_nolabels$`Cohort 2`,data_nolabels$`Cohort 3`))
y.train <- as.matrix(rbind(dist_labels$`Cohort 2`, dist_labels$`Cohort 3`))

x.test <- as.matrix(data_nolabels$`Cohort 1`)
y.testa <- as.matrix(dist_labels$`Cohort 1`)

#WE USED THE COMMENTED CODE BELOW TO DETERMINE THE BEST ALPHA VALUE - UNCOMMENT IF YOU 
#WISH TO RUN THIS FOR YOURSELF

## if using social distance, use gaussian settings
# fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
# fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)
# fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5) 
# fit.lasso.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1,
#                           family="gaussian")
# fit.ridge.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0,
#                           family="gaussian")
# fit.elnet.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
#                           family="gaussian")
# 
# for (i in 0:10) {
#   assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse",
#                                             alpha=i/10,family="gaussian"))
# }
for (i in 0) {
  assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse",
                                            alpha=i/10,family="gaussian"))
}

# ## Plot solution path and cross-validated MSE as function of lambda
# yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
# yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test)
# yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test)
# yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test)
# yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test)
# yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)
# yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test)
# yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test)
# yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test)
# yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test)
# yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
# 
# 
# #the following tests will only be relevant for gaussian models
# mse0 <- mean((y.testa - yhat0)^2)
# mse1 <- mean((y.testa - yhat1)^2)
# mse2 <- mean((y.testa - yhat2)^2)
# mse3 <- mean((y.testa - yhat3)^2)
# mse4 <- mean((y.testa - yhat4)^2)
# mse5 <- mean((y.testa - yhat5)^2)
# mse6 <- mean((y.testa - yhat6)^2)
# mse7 <- mean((y.testa- yhat7)^2)
# mse8 <- mean((y.testa - yhat8)^2)
# mse9 <- mean((y.testa - yhat9)^2)
# mse10 <- mean((y.testa - yhat10)^2)

#for gaussian data model prediction...
yhat0a <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
mse0_a <- mean((y.testa - yhat0a)^2)

library(ggpubr)
df <- as.data.frame(cbind(yhat0a,y.testa));names(df) <- c("pred","obs"); n <- length(df$pred)
dfa <- df
mod1plot <- ggscatter(df,x = "pred", y = "obs", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 3, label.y=3.5, label.sep = "\n"),
          xlab="predicted - ridge",ylab="", color = "darkslategrey",ggtheme = theme_pubclean(),ylim = c(1, 4),xlim=c(1,4))

R <- cor(yhat0a,y.testa,method="pearson")
R2 <- 1-(sum(unlist(lapply(1:n,function(i){(y.testa[i]-yhat0a[i])^2})))/sum(unlist(lapply(1:n,function(i){(y.testa[i]-mean(y.testa))^2}))))
print(paste("R = ",R,"R2 = ",R2))


########################
#model b
x.train <- as.matrix(rbind(data_nolabels$`Cohort 1`,data_nolabels$`Cohort 3`))
y.train <- as.matrix(rbind(dist_labels$`Cohort 1`, dist_labels$`Cohort 3`))

x.test <- as.matrix(data_nolabels$`Cohort 2`)
y.testb <- as.matrix(dist_labels$`Cohort 2`)

## if using social distance, use gaussian settings
# fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
# fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)
# fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5) 
# fit.lasso.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1,
#                           family="gaussian")
# fit.ridge.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0,
#                           family="gaussian")
# fit.elnet.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
#                           family="gaussian")
# 
# for (i in 0:10) {
#   assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse",
#                                             alpha=i/10,family="gaussian"))
# }
for (i in 0) {
  assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse",
                                            alpha=i/10,family="gaussian"))
}
# yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
# yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test)
# yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test)
# yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test)
# yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test)
# yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)
# yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test)
# yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test)
# yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test)
# yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test)
# yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
# 
# 
# #the following tests will only be relevant for gaussian models
# mse0 <- mean((y.testb - yhat0)^2)
# mse1 <- mean((y.testb - yhat1)^2)
# mse2 <- mean((y.testb - yhat2)^2)
# mse3 <- mean((y.testb - yhat3)^2)
# mse4 <- mean((y.testb - yhat4)^2)
# mse5 <- mean((y.testb - yhat5)^2)
# mse6 <- mean((y.testb - yhat6)^2)
# mse7 <- mean((y.testb- yhat7)^2)
# mse8 <- mean((y.testb - yhat8)^2)
# mse9 <- mean((y.testb - yhat9)^2)
# mse10 <- mean((y.testb - yhat10)^2)
## Plot solution path and cross-validated MSE as function of lambda

#for gaussian data model prediction...
yhat0b <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
mse0_b <- mean((y.testb - yhat0b)^2)

df <- as.data.frame(cbind(yhat0b,y.testb));names(df) <- c("pred","obs"); n <- length(df$pred)
dfb <- df
mod2plot <-ggscatter(df,x = "pred", y = "obs", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 3, label.y=3.5, label.sep = "\n"),
          xlab="predicted - ridge",ylab="observed", color=c("darkslategrey"),ggtheme = theme_pubclean(),ylim = c(1, 4),xlim=c(1,4))
R <- cor(yhat0b,y.testb,method="pearson")
R2 <- 1-(sum(unlist(lapply(1:n,function(i){(y.testb[i]-yhat0b[i])^2})))/sum(unlist(lapply(1:n,function(i){(y.testb[i]-mean(y.testb))^2}))))
print(paste("R = ",R,"R2 = ",R2))

########################
#model c
x.train <- as.matrix(rbind(data_nolabels$`Cohort 1`,data_nolabels$`Cohort 2`))
y.train <- as.matrix(rbind(dist_labels$`Cohort 1`, dist_labels$`Cohort 2`))

x.test <- as.matrix(data_nolabels$`Cohort 3`)
y.testc <- as.matrix(dist_labels$`Cohort 3`)

## if using social distance, use gaussian settings
# fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1)
# fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0)
# fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5) 
# fit.lasso.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=1,
#                           family="gaussian")
# fit.ridge.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=0,
#                           family="gaussian")
# fit.elnet.cv <- cv.glmnet(x.train, y.train, type.measure="mse", alpha=.5,
#                           family="gaussian")
# 
# for (i in 0:10) {
#   assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse",
#                                             alpha=i/10,family="gaussian"))
# }
for (i in 10) {
  assign(paste("fit", i, sep=""), cv.glmnet(x.train, y.train, type.measure="mse",
                                            alpha=i/10,family="gaussian"))
}

## Plot solution path and cross-validated MSE as function of lambda
# yhat0 <- predict(fit0, s=fit0$lambda.1se, newx=x.test)
# yhat1 <- predict(fit1, s=fit1$lambda.1se, newx=x.test)
# yhat2 <- predict(fit2, s=fit2$lambda.1se, newx=x.test)
# yhat3 <- predict(fit3, s=fit3$lambda.1se, newx=x.test)
# yhat4 <- predict(fit4, s=fit4$lambda.1se, newx=x.test)
# yhat5 <- predict(fit5, s=fit5$lambda.1se, newx=x.test)
# yhat6 <- predict(fit6, s=fit6$lambda.1se, newx=x.test)
# yhat7 <- predict(fit7, s=fit7$lambda.1se, newx=x.test)
# yhat8 <- predict(fit8, s=fit8$lambda.1se, newx=x.test)
# yhat9 <- predict(fit9, s=fit9$lambda.1se, newx=x.test)
# yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
# 
# 
# #the following tests will only be relevant for gaussian models
# mse0 <- mean((y.testc - yhat0)^2)
# mse1 <- mean((y.testc - yhat1)^2)
# mse2 <- mean((y.testc - yhat2)^2)
# mse3 <- mean((y.testc - yhat3)^2)
# mse4 <- mean((y.testc - yhat4)^2)
# mse5 <- mean((y.testc - yhat5)^2)
# mse6 <- mean((y.testc - yhat6)^2)
# mse7 <- mean((y.testc- yhat7)^2)
# mse8 <- mean((y.testc - yhat8)^2)
# mse9 <- mean((y.testc - yhat9)^2)
# mse10 <- mean((y.testc - yhat10)^2)
#for gaussian data model prediction...
yhat1c <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
mse10_c <- mean((y.testc - yhat1c)^2)


df <- as.data.frame(cbind(yhat1c,y.testc));names(df) <- c("pred","obs"); n <- length(df$pred)
dfc <- df
mod3plot <- ggscatter(df,x = "pred", y = "obs", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.coeff.args = list(method = "pearson", label.x = 3, label.y=3.5,label.sep = "\n"),
          xlab="predicted - lasso",ylab="", color=c("darkslategrey"),ggtheme = theme_pubclean(),ylim = c(1, 4),xlim=c(1,4))
R <- cor(yhat1c,y.testc,method="pearson")
R2 <- 1-(sum(unlist(lapply(1:n,function(i){(y.testc[i]-yhat1c[i])^2})))/sum(unlist(lapply(1:n,function(i){(y.testc[i]-mean(y.testc))^2}))))
print(paste("R = ",R,"R2 = ",R2))

#####################################
#save dfa, dfb and dfc for use in plotting in case you want to do it later
save(dfa,dfb,dfc , file = "dfs.RData")

#####################################
#open the df files created previously
library(ggplot2)
library(ggpubr)
load("dfs.RData")
df_list <- list(dfa, dfb, dfc)
names(df_list) <- c('dfa', 'dfb', 'dfc')
# #concatenate
# all_yhats <- rbind(yhat0a, yhat0b, yhat1c)
# all_ytests <- rbind(y.testa,y.testb,y.testc)
# 
# df <- as.data.frame(cbind(all_yhats,all_ytests));names(df) <- c("pred","obs"); n <- length(df$pred)
# ggscatter(df,x = "pred", y = "obs", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson",
#           xlab="predicted class - alpha=0",ylab="observed class")
# R <- cor(all_yhats,all_ytests,method="pearson")
# R2 <- 1-(sum(unlist(lapply(1:n,function(i){(all_ytests[i]-all_yhats[i])^2})))/sum(unlist(lapply(1:n,function(i){(all_ytests[i]-mean(all_ytests))^2}))))
# print(paste("R = ",R,"R2 = ",R2))
################
#meta-analysis


data_a <- as.data.frame(cbind(y.testa,yhat0a,net_corr_vectors$`Cohort 1`$wholebrain$ppt_1,net_corr_vectors$`Cohort 1`$wholebrain$ppt_2))
colnames(data_a) <- c("y.test","yhat","ppt_1","ppt_2")
model_a <- lmer(y.test ~ yhat +(1|ppt_1)+(1|ppt_2), data_a)#y.test is dependent variable, yhat independent variable
# model_a <- lmer(y.test ~ yhat +(1+yhat||ppt_1)+(1+yhat||ppt_2), data_a)#with random slopes

data_b <- as.data.frame(cbind(y.testb,yhat0b,net_corr_vectors$`Cohort 2`$wholebrain$ppt_1,net_corr_vectors$`Cohort 2`$wholebrain$ppt_2))
colnames(data_b) <- c("y.test","yhat","ppt_1","ppt_2")
model_b <- lmer(y.test ~ yhat +(1|ppt_1)+(1|ppt_2), data_b)
# model_b <- lmer(y.test ~ yhat +(1+yhat||ppt_1)+(1+yhat||ppt_2), data_b)#with random slopes

data_c <- as.data.frame(cbind(y.testc,yhat1c,net_corr_vectors$`Cohort 3`$wholebrain$ppt_1,net_corr_vectors$`Cohort 3`$wholebrain$ppt_2))
colnames(data_c) <- c("y.test","yhat","ppt_1","ppt_2")
model_c <- lmer(y.test ~ yhat +(1|ppt_1)+(1|ppt_2), data_c)
# model_c <- lmer(y.test ~ yhat +(1+yhat||ppt_1)+(1+yhat||ppt_2), data_c)#with random slopes

model_list <- list(model_a,model_b,model_c)
names(model_list) <- c("model 1", "model 2", "model 3")

MetaAnalysis <- meta_models(model_list = model_list)
metaforlmer::ggforest(MetaAnalysis, labels = "predictive model performance")
#######################################
#create a list with the models in it and get a p value for each of the model t values
mdls <- list(model_a,model_b,model_c)
names(mdls) <- c("model_a","model_b","model_c")

pvals <- list()
for (mdl in names(mdls)){
  pv <- summary(as_lmerModLmerTest(model = mdls[[mdl]]))
  pvals[[mdl]] <- pv$coefficients["yhat","Pr(>|t|)"]
}
# pval_a <- summary(as_lmerModLmerTest(model = model_a))
# p_a <- pval_a$coefficients["yhat","Pr(>|t|)"]
# pval_b <- summary(as_lmerModLmerTest(model = model_b))
# p_b <- pval_b$coefficients["yhat","Pr(>|t|)"]
# pval_c <- summary(as_lmerModLmerTest(model = model_c))
# p_c <- pval_c$coefficients["yhat","Pr(>|t|)"]
# 


cors <- list()
for (d in names(df_list)){
  cors[[d]] <- cor(df_list[[d]]$pred,df_list[[d]]$obs,method = c("pearson"))
}
#now you have two lists: pvals and cors. These contain the correlations and pvalues you will add to the plot you create next.
dfa$model <- c("model 1");dfb$model <- c("model 2"); dfc$model <- c("model 3")
df_all <- rbind(dfa,dfb,dfc)
names(df_all) <- c("predicted","observed", "model")

r1 <- c(paste("Model 1: R =",gsub("0\\.",".",round(cors[["dfa"]], digits = 3)),"; p =",gsub("0\\.",".",round(pvals[["model_a"]], digits = 3))))
r2 <- c(paste("Model 2: R =",gsub("0\\.",".",round(cors[["dfb"]], digits = 3)),"; p =",gsub("0\\.",".",round(pvals[["model_b"]], digits = 3))))
r3 <- c(paste("Model 3: R =",gsub("0\\.",".",round(cors[["dfc"]], digits = 3)),"; p =",gsub("0\\.",".",round(pvals[["model_c"]], digits = 3))))

#create a matrix where x=y

#library(viridis)
pal <- c("#0b409c","#ff0000","#fd5f00","#f8b500")
#pal <- magma(4)
# pal <- c('mediumaquamarine',"darkpink",'darkgoldenrod1', "darkslategrey")
p <- ggplot(data=df_all, aes(x=observed, y=predicted, color=model)) + 
  geom_point(alpha = .5,size=1.5, position = position_dodge(width=0.1)) + geom_smooth(method=lm, aes(fill=model))+
  xlab(c("Observed distance"))+ ylab(c("Predicted distance"))+
  scale_color_manual("Prediction performance",values = pal, labels=c(r1,r2,r3)) + 
  scale_fill_manual("Prediction performance",values = pal, labels=c(r1,r2,r3))+
  theme_pubclean()+ylim(1,4)
  p <- p + theme(legend.position = c(0.33,0.8),
                 axis.ticks.y = element_blank())+ 
    geom_abline(intercept = 0, slope = 1,color="lightblue", 
                linetype="dashed", size=1)



#####################

names(model_list) <- c("Model 1", "Model 2", "Model 3") 
source(paste(home_dir,'ggforest_try_col.R',sep = ""))
MetaAnalysis <- meta_models(model_list = model_list)
plot_meta <- ggforest(MetaAnalysis, labels = "Predictive model performance", palette = pal)
ggarrange(p,plot_meta, ncol = 2, labels = c("a", "b"), common.legend = F,
          nrow = 1) 


#sqrt(mse0_a);sqrt(mse0_b);sqrt(mse10_c)
tiff(paste(home_dir,c("fig8.tiff"),sep = ""), width = 2126, height = 1400, units = "px", res = 300)
# 2. Create the plot

ggarrange(p,plot_meta, ncol = 2, labels = c("a", "b"), common.legend = F,
          nrow = 1) 
# 3. Close the file
dev.off()
