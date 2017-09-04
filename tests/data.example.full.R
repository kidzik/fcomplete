library("devtools") ; library("roxygen2") ; library("fpca")
roxygenise() ; install(".") ; library("fcomplete")
# if (!("all.data" %in% ls()))
#   all.data = read.csv("/media/lukasz/SECURED/annotations/data.csv")
if (!("all.data" %in% ls())){
#  all.data = read.csv("/home/lukasz/Dropbox/DATA/CP/alldata.csv")
#  gait.cycles = t(read.csv("/home/lukasz/Dropbox/DATA/CP/G_avg_CP.csv"))
  all.data = read.csv("/home/lukasz/alldata.csv")
  gait.cycles = t(read.csv("/home/lukasz/G_avg_CP.csv"))
}
pcas = prcomp(gait.cycles)
all.data = cbind(all.data, pcas$x[,1:10])

# all.data[all.data == "NaN"] = NA
# all.data = all.data[order(all.data[,c("age")]),]
# all.data = all.data[order(all.data[,c("Patient_ID")]),]
all.data.subset = all.data[all.data$side == "L",]
all.data.subset = all.data.subset[all.data.subset$age < 15,] # remove outliers
# all.data = all.data[all.data$side == "L",]
# all.data = all.data[(all.data$age) < 18,] # remove outliers

# all.data$mass[is.na(all.data$mass)] = mean(all.data$mass,na.rm=TRUE)
# all.data$height[is.na(all.data$height)] = mean(all.data$height,na.rm=TRUE)
# all.data$cadence[is.na(all.data$cadence)] = mean(all.data$cadence,na.rm=TRUE)
# all.data$bmi = all.data$mass / all.data$height**2
cmeans = colMeans(all.data.subset[,300:309])
all.data.subset$score = sqrt(rowSums(t(t(all.data.subset[,300:309]) - cmeans)**2))

all.data.filtered = all.data.subset[,c("Patient_ID", "speed","cadence","age","bmi","height","KneeFlex_maxExtension","O2cost","PC1","PC2","PC3")]
all.data.filtered = all.data.filtered[order(all.data.filtered$Patient_ID, all.data.filtered$age),]
TBL = table(all.data.filtered$Patient_ID)
all.data.filtered = all.data.filtered[all.data.filtered$Patient_ID %in% names(TBL[TBL > 2]),]
all.data.filtered = na.omit(all.data.filtered)
all.data.filtered = all.data.filtered[!is.nan(all.data.filtered$PC1),]

data = sample.long(all.data.filtered, "Patient_ID", "age", "PC1", ratio = 0.05)

# IMPUTE
#model.impute = fregression(GDI_Pre:age ~ 1 | Patient_ID, data$train, lambda= c(100), thresh = 1e-6)
#model.impute = fregression(GDI_Pre:age ~ 1 | Patient_ID, data$train, lambda= c(0,1,5,10,25,50,100,500,1000,5000), thresh = 1e-4,method = "fimpute", K=3)
roxygenise() ; install(".") ; library("fcomplete")

lambdas = seq(0,10,length.out = 20)
lambdas.reg = seq(0,10,length.out = 20)
model.regression = fregression(PC1:age ~ PC1 + bmi + O2cost | Patient_ID, data$train, method = "fimpute", K = 5, thresh=1e-5, K.reg = 5, lambda = lambdas, lambda.reg = lambdas.reg, d=6)
sqrt(mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE)) #/ ss

model.regression1 = fregression(PC1:age ~ speed + bmi + PC1 | Patient_ID, data$train, method = "fimpute", K = 6, thresh=1e-5, K.reg = 18, lambda = lambdas, lambda.reg = lambdas.reg, d=6)
model.regression = fregression(PC1:age ~ PC1 | Patient_ID, data$train, method = "fimpute", K = 6, thresh=1e-5, K.reg = 6, lambda = lambdas, lambda.reg = lambdas.reg, d=6)
sqrt(mean((model.regression1$fit - data$test.matrix)**2, na.rm = TRUE)) #/ ss
sqrt(mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE)) #/ ss

#model.impute = fregression(PC1:age ~ 1 | Patient_ID, data$train, lambda= c(1.5), thresh = 1e-6, method = "fimpute", K=2, d=6)
model.impute.fpcs = fregression(PC1:age ~ 1 | Patient_ID, data$train, lambda= c(7.5), thresh = 1e-4, method = "fpcs", K=2:6, d=6)
model.mean = fregression(PC1:age ~ 1 | Patient_ID, data$train, lambda= c(100), thresh = 1e-6, method = "mean", K=2)

# REGRESSION
# model.regression = fregression(GDI_Pre:age ~ speed + mass + height + cadence| Patient_ID, all.data, method = "fpcs", K = 3,thresh=1e-4)

 #res = test.long(data, model.impute$fit, as.numeric(colnames(model.regression$Y)), rownames(model.regression$Y))
#cor(res)

# DIM REDUCTION
#model.unsupervised = fregression(1:age ~ speed + leglen + height | Patient_ID, all.data)

ss = sqrt(mean((data$test.matrix)**2, na.rm = TRUE))

#lambdas = (0:20)/4
lambdas = seq(2,7,length.out = 20) # (0 + (0:10)*1)
model.impute = fregression(PC1:age ~ 1 | Patient_ID, data$train, lambda= lambdas, thresh = 1e-5, method = "fimpute", K=6, d=6)
model.impute$meta

sqrt(mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE)) #/ ss
sqrt(mean((model.impute$fit - data$test.matrix)**2, na.rm = TRUE)) #/ ss
sqrt(mean((model.impute.fpcs$fit - data$test.matrix)**2, na.rm = TRUE)) #/ ss
sqrt(mean((model.mean$fit - data$test.matrix)**2, na.rm = TRUE)) #/ ss

ind = row.names(data$test.matrix) %in% data$X$Patient_ID[data$test.ob[1:3]]

#matplot(t(model.regression$Y[ind,]),t='l',lty=1,lwd=4,ylim = lims)
matplot(t(data$test.matrix[ind,]),t='p',lty=2,lwd=2,pch="x", ylim = c(-200,200))
matplot(t(model.regression$fit[ind,]),t='l',lty=4,add=T,lwd=2)
matplot(t(model.mean$fit[ind,]),t='l',lty=3,add=T,lwd=2)
matplot(t(model.impute$fit[ind,]),t='l',lty=2,add=T,lwd=2)
matplot(t(model.impute.fpcs$fit[ind,]),t='l',lty=1,add=T,lwd=2)
title("Mean prediction")
