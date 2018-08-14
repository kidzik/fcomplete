# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools")
devtools::install(".")
library("fcomplete")
library("ggplot2")
library("latex2exp")

res = list()
nexp = 1
dgrid = 31
d = 7

# SIMULATE DATA
test.experiment = function(exp.id){
set.seed(30 + exp.id)
simulation = fsimulate(dgrid = dgrid,clear = 0.9, n = 100, noise.mag = 0.3, d = d, K = 2)
data = simulation$data
ftrue = simulation$ftrue
K = simulation$params$K

# TUNING PARAMS
lambdas.pca = seq(0,2,length.out = 5)
lambdas.reg = seq(0,1,length.out = 5)

# model.fslr = fregression(Y:time ~ Y + X1 + X2 | id, data, lambda = lambdas.pca, thresh = 1e-4, lambda.reg = lambdas.reg, method = "fimpute", bins = dgrid)

# model.mean = fregression(Y:time ~ 1 | id, data, method = "mean", bins = dgrid)
model.fpca = fregression(Y:time ~ 1 | id, data, lambda = 0, K = 2:d, thresh = 1e-7, method = "fpcs", bins = dgrid)
model.fimp = fregression(Y:time ~ 1 | id, data, lambda = lambdas.pca, thresh = 0, final = "soft", maxIter = 1000, fold = 5, cv.ratio = 0.05, bins = dgrid)

# REPORT RESULTS
errors = c(
  mean((ftrue - mean(data$Y))**2),
  mean((ftrue - model.fpca$fit)**2),
  mean((ftrue - model.fimp$fit)**2)
)
}
res = lapply(1:1, test.experiment)

errors = c()
for (i in 1:length(res)){
  errors = rbind(errors, res[[i]])
}
print(errors)
colMeans(errors)

tbl.true = cbind(
  errors,
  100*(1-errors/errors[1])
)
colnames(tbl.true) = c("MSE","% expl")
rownames(tbl.true) = c("mean",
                       "fpca",
                       "fimpute",
                       "regression")
print(tbl.true)

# SAVE EXPERIMENT RESULTS
res = list()
res$tbl = tbl.true
res$errors = errors
res$mean = model.mean
res$fimp = model.fimp
res$fpca = model.fpca
res$fslr = model.fslr
res$simulation = simulation
res
