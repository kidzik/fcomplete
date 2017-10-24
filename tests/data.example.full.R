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

res = c()
res = read.csv("res.csv",row.names = 1)
par(cex=1.3)
boxplot(t(res))

for (i in 1:10){
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

errors = c(sqrt(mean((model.regression$fit - data$test.matrix)**2, na.rm = TRUE)),
  sqrt(mean((model.impute$fit - data$test.matrix)**2, na.rm = TRUE)),
  sqrt(mean((model.impute.fpcs$fit - data$test.matrix)**2, na.rm = TRUE)),
  sqrt(mean((model.mean$fit - data$test.matrix)**2, na.rm = TRUE)))
names(errors) = c("regression","impute","fPCA","mean")

res = cbind(res, errors)
}
rownames(res) = c("regression","impute","fPCA","mean")
cbind(rowMeans(res),
apply(res,FUN=sd,1))
colnames(res) = paste("run",1:11)
par(mfrow=c(1,1))
boxplot(t(res[c(3,2,1,4),]))

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
library(RColorBrewer)

#dd = all.data.filtered[,c("Patient_ID","age","PC1")]

dd = data$train[,c("Patient_ID","age","PC1")] #fcomplete:::fc.wide2long(model.impute$data[[1]]$train)
#colnames(dd) = c("Patient_ID","age","PC1")
dd$Patient_ID = as.factor(dd$Patient_ID)

pp = ggplot(aes(x = age, y = PC1, color = Patient_ID), data = dd[1:200,]) + ylab("1st component") +
  geom_point(size = 3) + theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  stat_function(fun = approxfun(lowess(dd$age,dd$PC1)), size = 2.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
pp
ggsave("~/Dropbox/Presentations/Mobilize17/images/fcomplete/points.pdf")
pp + geom_line(size=0.7)
ggsave("~/Dropbox/Presentations/Mobilize17/images/fcomplete/grouped.pdf")

# The palette with grey:
# The palette with black:
cbPalette <- c("#000000", "#A69F00", "#56B4E9", "#009E73", "#008442") #, "#0072B2", "#D55E00", "#CC79A7")

gg = as.numeric(colnames(model.impute.fpcs$Y))
pp = ggplot(aes(x = age, y = PC1, color = Patient_ID), data = dd[2:14,]) +
  scale_fill_manual(values = cbPalette) + scale_colour_manual(values = cbPalette) +
  ylab("1st component") +
  geom_point(size = 3) + theme_set(theme_grey(base_size = 26)) + theme(legend.position="none", panel.background = element_rect(fill = "white",linetype = 1,colour = "grey50",size = 1,)) +
  stat_function(fun = approxfun(lowess(dd$age,dd$PC1)), size = 2.5, colour = "#000000")+ scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) + xlim(5,13)
pp + geom_line(size = 0.7) +
  stat_function(fun = approxfun(gg, model.impute.fpcs$fit[5,]) , color = cbPalette[5]) +
  stat_function(fun = approxfun(gg, model.impute.fpcs$fit[2,]) , color = cbPalette[2]) +
  stat_function(fun = approxfun(gg, model.impute.fpcs$fit[3,]) , color = cbPalette[3])  +
  stat_function(fun = approxfun(gg, model.impute.fpcs$fit[4,]) , color = cbPalette[4])

#pp #+ model.impute.fpcs$fit[1:4,]



dgrid = 51
d = 6
basis = fda::create.bspline.basis(c(0,1), d)
plot(basis)
S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis)

par(mfrow=c(1,2),cex=1.3)
matplot(0:50/50, (S),t='l',ylab = "value",xlab="time",ylim = c(-1,1), lwd=4, col = brewer.pal(d, "Paired"), lty=1)
title("Spline basis")
matplot(0:50/50, (svd(S)$u),t='l',ylab="value",xlab="time",ylim = c(-1,1), lwd=4, col = brewer.pal(d, "Paired"),lty=1)
title("Orthonormalized spline basis")
