# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools") ; devtools::install(".")
library("fcomplete")

#rm(list = ls())
res = list()
nexp = 1
dgrid = 51

for(exp.id in 1:nexp){

# SIMULATE DATA
set.seed(323 + exp.id)
simulation = fsimulate(dgrid = dgrid,clear = 0.85)
data = simulation$data
ftrue = simulation$ftrue
K = 6 #simulation$params$K

# REGRESSION
model.mean = fregression(Y:time ~ 1 | id, data, method = "mean", bins = dgrid)
model.fpca = fregression(Y:time ~ 1 | id, data, lambda = 0, K = 2:K, thresh = 1e-7, method = "fpcs", bins = dgrid)

lambdas = seq(0,10,length.out = 20)
model.fimp = fregression(Y:time ~ 1 | id, data, lambda = lambdas, thresh = 0, final = "soft", maxIter = 100, fold = 5, cv.ratio = 0.05, K = K, bins = dgrid)
model.fcmp = fregression(0:time ~ Y + X1 + X2 | id, data, lambda = lambdas, K = K, final = "soft")
lambdas = seq(0,2,length.out = 10)
model.freg = fregression(Y:time ~ X1 + X2 | id, data, lambda = lambdas, thresh = 1e-5, lambda.reg = 0.1 * 5:20, method = "fimpute", K = K)

# REPORT RESULTS
idx = unique(data$id)
errors = c(
  mean((ftrue[idx,] - model.mean$fit)**2),
  mean((ftrue[idx,] - model.fpca$fit)**2),
  mean((ftrue[idx,] - model.fimp$fit)**2),
  mean((ftrue[idx,] - model.fcmp$fit)**2),
  mean((ftrue[idx,] - model.freg$fit)**2)
)
errors
model.fimp$meta

tbl.true = cbind(
  errors,
  100*(1-errors/errors[1])
)
colnames(tbl.true) = c("MSE","% expl")
rownames(tbl.true) = c("mean","fpca","fimpute","fcompress","regression")
print(tbl.true)

res[[exp.id]] = list()
res[[exp.id]]$tbl = tbl.true
res[[exp.id]]$fimpute = model.fimp
res[[exp.id]]$fpca = model.fpca
res[[exp.id]]$freg = model.fimp

# PLOT EXAMPLES
par(mfrow=c(2,2), cex=1.3)
ind = 1:2 + 10
idx = as.numeric(model.freg$id)

lims = c(min(ftrue[idx,][ind,],simulation$fobs[ind,],na.rm = TRUE) - 2,
         max(ftrue[idx,][ind,],simulation$fobs[ind,],na.rm = TRUE) + 2)
matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.mean$fit[ind,]),t='l',lty=2,add=T,lwd=2)
title("Mean prediction")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.fimp$fit[ind,]),t='l',lty=2,add=T,lwd=2)
title("Functional impute")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.fpca$fit[ind,]),t='l',lty=2, add=T,lwd=2)
title("Functional PCA")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims,xlab = "time")
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.freg$fit[ind,]),t='l',lty=2, add=T,lwd=2)
title("Functional regression")

par(mfrow=c(1,2), cex=1.3)
matplot(t(-model.fimp$v)[,1:3],t='l',lwd = 4)
title("First 3 PCs from fPCA method")
matplot(t(model.fpca$v)[,1:3],t='l', lwd = 4)
title("First 3 singular vectors from SFI method")

par(mfrow=c(1,2), cex=1.3)
plot(cv.err ~ lambda, model.fimp$meta, type='o', ylab="CV error")
title("Sample error by cross-validation")
model.fimp$meta$cv.K = round(model.fimp$meta$cv.K)
plot(cv.K ~ lambda, model.fimp$meta, ylab="CV rank")
title("Rank by cross-validation")

n = nrow(model.fcmp$u)
#plot(model.fcmp$u[,1:2], col = 1 + (1:n > n/3), xlab = "Component 1", ylab="component 2")
#title("Decomposition of (Y(t),X_1(t),X_2(t))")
}

joint.tbl = res[[1]]$tbl
for (i in 2:length(res)){
  joint.tbl = joint.tbl + res[[i]]$tbl
}
joint.tbl = joint.tbl/length(res)
res[[1]]$fimpute$meta


# PLOT model.fimp
par(mfrow=c(1,1))
plot(model.fimp$meta[1:100,c(1,2)],ylim=c(0,10))
lines(model.fimp$meta[1:100,c(1,3)] )
