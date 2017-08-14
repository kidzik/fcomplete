rm(list = ls())

# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools") ; devtools::install(".")
library("fcomplete")

# SIMULATE DATA
set.seed(1)
simulation = fsimulate()
data = simulation$data
ftrue = simulation$ftrue
K = simulation$params$K

# REGRESSION
model.mean = fregression(Y:time ~ 1 | id, data, method = "mean")
model.fpca = fregression(Y:time ~ 1 | id, data, lambda = 0, K = c(3,4,5), thresh = 1e-7, method = "fpcs")

lambdas = c(2,3,4,5,6,8,10,12,15,20)
model.fimp = fregression(Y:time ~ 1 | id, data, lambda = lambdas, thresh = 1e-5, final = "hard")
model.fcmp = fregression(0:time ~ Y + X1 + X2 | id, data, lambda = lambdas, K = K, final = "hard")
model.freg = fregression(Y:time ~ X1 + X2 | id, data, lambda = lambdas, thresh = 1e-5, lambda.reg = 0.1, method = "fpcs", K = K)

# REPORT RESULTS
idx = unique(data$id)
errors = c(
  sqrt(mean((ftrue[idx,] - model.mean$fit)**2)),
  sqrt(mean((ftrue[idx,] - model.fpca$fit)**2)),
  sqrt(mean((ftrue[idx,] - model.fimp$fit)**2)),
  sqrt(mean((ftrue[idx,] - model.fcmp$fit[[1]])**2)),
  sqrt(mean((ftrue[idx,] - model.freg$fit)**2))
)

tbl.true = cbind(
  errors,
  100*(1-errors/errors[1])
)
colnames(tbl.true) = c("MSE","% expl")
rownames(tbl.true) = c("mean","fpca","fimpute","fcompress","regression")
print(tbl.true)

# PLOT EXAMPLES
par(mfrow=c(2,2))
ind = 10 + 1:2
idx = as.numeric(model.freg$id)

lims = c(min(simulation$fobs[ind,],na.rm = TRUE) - 0.5,
         max(simulation$fobs[ind,],na.rm = TRUE) + 0.5)
matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims)
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.mean$fit[ind,]),t='l',lty=2,add=T,lwd=2)
title("Mean prediction")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims)
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.fimp$fit[ind,]),t='l',lty=2,add=T,lwd=2)
title("Functional impute")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims)
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.fpca$fit[ind,]),t='l',lty=2, add=T,lwd=2)
title("Functional SVD")

matplot(t(ftrue[idx,][ind,]),t='l',lty=1,lwd=4,ylim = lims)
matplot(t(simulation$fobs[ind,]),t='p',lty=2,add=T,lwd=2,pch="x")
matplot(t(model.freg$fit[ind,]),t='l',lty=2, add=T,lwd=2)
title("Functional regression")

par(mfrow=c(1,1), cex=1.3)
plot(cv.err ~ lambda, model.fimp$meta, type='o', ylab="CV error")
title("Cross-validation of sample error")
plot(cv.K ~ lambda, model.fimp$meta, ylab="CV dimensionality")
title("Cross-validation of the dimension")

n = nrow(model.fcmp$u)
plot(model.fcmp$u[,1:2], col = 1 + (1:n > n/3), xlab = "Component 1", ylab="component 2")
title("Decomposition of (Y(t),X_1(t),X_2(t))")
