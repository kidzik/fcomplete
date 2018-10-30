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
lambdas.reg = 0
lambdas.pca = 1

# model.mean = fregression(Y:time ~ 1 | id, data, method = "mean", bins = dgrid)
model.fpca = fregression(Y:time ~ 1 | id, data, lambda = 0, K = 2:d, thresh = 1e-7, method = "fpcs", bins = dgrid, maxIter = 1000)
model.fimp = fregression(Y:time ~ 1 | id, data, lambda = lambdas.pca, thresh = 0, final = "soft", maxIter = 1000, fold = 5, cv.ratio = 0.05, bins = dgrid)

model.fslr = fregression(Y:time ~ Y + X1| id, data, K=2, K.reg=2, lambda = lambdas.pca, thresh = 1e-10, lambda.reg = lambdas.reg, method = "fimpute", bins = dgrid, maxIter = 5000)

# REPORT RESULTS
errors = c(
  mean((ftrue - mean(data$Y))**2),
  mean((ftrue - model.fpca$fit)**2),
  mean((ftrue - model.fimp$fit)**2),
  mean((ftrue - model.fslr$fit)**2)
)
errors
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
