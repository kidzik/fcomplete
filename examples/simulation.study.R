# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools")
devtools::install(".")
library("fcomplete")
library("ggplot2")
library("latex2exp")
library("parallel")

res = list()
nexp = 1
dgrid = 51
d = 7

# SIMULATE DATA
#test.experiment = function(exp.id){
set.seed(30 + exp.id)
#simulation = fsimulate(dgrid = dgrid,clear = 0.9, n = 50, noise.mag = 0.3, d = d, K = 2)
simulation = fsimulate(dgrid = dgrid,num_points = 3, n = 100, noise.mag = 0.3, d = d, K = 2)
data = simulation$data
ftrue = simulation$ftrue
K = simulation$params$K

# TUNING PARAMS
lambdas.pca = seq(1.5,3,length.out = 10) / (dgrid)
lambdas.reg = seq(0,0.2,length.out = 20)
#lambdas.reg = 0
#lambdas.pca = 1

model.mean = fregression(Y ~ time | id, data, method = "mean", bins = dgrid)
#model.fpca = fregression(Y ~ time | id, data, lambda = 0, K = 2:d, thresh = 1e-7, method = "fpcs", bins = dgrid, maxIter = 1000)
model.fimp = fregression(Y ~ time | id, data, lambda = lambdas.pca, thresh = 1e-2, final = "soft", maxIter = 1000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
mean(((ftrue - model.fimp$fit)[!is.na(simulation$fobs)])**2)

#model.X1 = fregression(X1 + X2 + Y ~ time | id, data, K=4, lambda = lambdas.pca, thresh = 1e-10, method = "fimpute", bins = dgrid, maxIter = 5000)
#model.fslr = fregression(Y ~ time + U1 + U2 + U3 + U4 | id, data, model.X1$u, K=3, thresh = 1e-10, lambda = lambdas.reg, method = "fimpute", bins = dgrid, maxIter = 5000)

# REPORT RESULTS
errors = c(
  mean((ftrue - mean(data$Y))**2),
  mean((ftrue - model.fpca$fit)**2),
  mean((ftrue - model.fimp$fit)**2)
#  mean((ftrue - model.fslr$fit)**2)
)
errors
#}
res = mclapply(1:5, test.experiment, mc.cores = 5)

errors = c()
for (i in 1:length(res)){
  errors = rbind(errors, res[[i]])
}
print(errors)
colMeans(errors)

colnames(errors) = c("mean",
                       "fpca",
                       "fimpute" #,
#                       "regression"
                       )
print(errors)
