library("devtools") ; library("roxygen2") ; library("fpca")
roxygenise() ; install(".") ; library("fcomplete")

#############################
# PREPARE DATA
#############################
if (!("all.data" %in% ls())){
  all.data = read.csv("/home/lukasz/Dropbox/DATA/CP/alldata.csv")
  gait.cycles = t(read.csv("/home/lukasz/Dropbox/DATA/CP/G_avg_CP.csv"))
#  all.data = read.csv("/home/lukasz/alldata.csv")
#  gait.cycles = t(read.csv("/home/lukasz/G_avg_CP.csv"))
}
pcas = prcomp(gait.cycles)
all.data = cbind(all.data, pcas$x[,1:10])

all.data.subset = all.data[all.data$side == "L",]
all.data.subset = all.data.subset[all.data.subset$age < 15,] # remove outliers
all.data.subset = all.data.subset[all.data.subset$bmi > 5,] # remove outliers

cmeans = colMeans(all.data.subset[,300:309])
all.data.subset$score = sqrt(rowSums(t(t(all.data.subset[,300:309]) - cmeans)**2))

all.data.filtered = all.data.subset[,c("Patient_ID", "speed","cadence","age","bmi","height","KneeFlex_maxExtension","O2cost","PC1","PC2","PC3")]
all.data.filtered = all.data.filtered[order(all.data.filtered$Patient_ID, all.data.filtered$age),]
TBL = table(all.data.filtered$Patient_ID)
all.data.filtered = all.data.filtered[all.data.filtered$Patient_ID %in% names(TBL[TBL > 2]),]
all.data.filtered = na.omit(all.data.filtered)
all.data.filtered = all.data.filtered[!is.nan(all.data.filtered$PC1),]
#############################

# Sample data for testing
data = sample.long(all.data.filtered, "Patient_ID", "age", "PC1", ratio = 0.05)

# Set up parameters
lambdas = seq(0,10,length.out = 21)
K = 6

# IMPUTE
model.impute = fregression(PC1:age ~ 1 | Patient_ID, data$train, lambda= lambdas, thresh = 1e-5, method = "fimpute", K=K, d=K)
model.impute.fpcs = fregression(PC1:age ~ 1 | Patient_ID, data$train, lambda= c(7.5), thresh = 1e-4, method = "fpcs", K=2:K, d=K)
model.mean = fregression(PC1:age ~ 1 | Patient_ID, data$train,  method = "mean")

# REGRESSION
model.regression = fregression(PC1:age ~ PC1 + bmi + O2cost | Patient_ID, data$train, method = "fimpute", K = K, thresh=1e-5, K.reg = K, lambda = lambdas, lambda.reg = lambdas, d=K)

ss = sqrt(mean((data$test.matrix - mean(data$test.matrix, na.rm = TRUE))**2, na.rm = TRUE))

sqrt(mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE)) / ss
sqrt(mean((model.impute$fit - data$test.matrix)**2, na.rm = TRUE)) / ss
sqrt(mean((model.impute.fpcs$fit - data$test.matrix)**2, na.rm = TRUE)) / ss
sqrt(mean((model.mean$fit - data$test.matrix)**2, na.rm = TRUE)) / ss

ind = row.names(data$test.matrix) %in% data$X$Patient_ID[data$test.ob[1:3]]

matplot(t(data$test.matrix[ind,]),t='p',lty=2,lwd=2,pch="x", ylim = c(-100,100))
matplot(t(model.regression$fit[ind,]),t='l',lty=4,add=T,lwd=2)
matplot(t(model.mean$fit[ind,]),t='l',lty=3,add=T,lwd=2)
matplot(t(model.impute$fit[ind,]),t='l',lty=2,add=T,lwd=2)
matplot(t(model.impute.fpcs$fit[ind,]),t='l',lty=1,add=T,lwd=2)
title("Mean prediction")


#############################
# OTHER PLOTS FOR PAPER
#############################
library("ggplot2")
library("ggthemes")

dd = all.data.filtered[,c("Patient_ID","age","bmi")]
dd$Patient_ID = as.factor(dd$Patient_ID)

pp = ggplot(aes(x = age, y = bmi, color = Patient_ID), data = dd[1:200,]) +
  geom_point() + theme_set(theme_grey(base_size = 18)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  stat_function(fun = approxfun(lowess(dd$age,dd$bmi)), size = 2.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

pp
pp + geom_line()
