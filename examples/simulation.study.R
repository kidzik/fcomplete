# Regenerate and reinstall fimpute
library("roxygen2") ; roxygenize()
library("devtools")
#devtools::install(".")
library("fcomplete")
library("ggplot2")
library("latex2exp")
library("parallel")

res = list()
dgrid = 51
d = 7

# SIMULATE DATA
test.experiment = function(exp.id){
#set.seed(35 + exp.id)
#simulation = fsimulate(dgrid = dgrid,clear = 0.9, n = 50, noise.mag = 0.3, d = d, K = 2)
simulation = fsimulate(dgrid = dgrid,num_points = 3, n = 100, noise.mag = 0.1, d = d, K = 4)
data = simulation$data
ftrue = simulation$ftrue
K = simulation$params$K

# TUNING PARAMS
lambdas.pca = seq(0,0.1,length.out = 5)
#lambdas.pca = 0
lr = 0.5
#lambdas.reg = seq(0,0.2,length.out = 20)
#lambdas.reg = 0
#lambdas.pca = 1

model.mean = fregression(Y ~ time | id, data, method = "mean", bins = dgrid)

ptm <- proc.time()
model.fpca = fregression(Y ~ time | id, data, lambda = 0, K = 2:d, thresh = 1e-4, method = "fpcs", bins = dgrid, maxIter = 1000)
time.fpca = (proc.time() - ptm)[3]

ptm <- proc.time()
model.fimp = fregression(Y ~ time | id, data, method="proximal_grad_grid", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
#model.fimp = fregression(Y ~ time | id, data, lambda = lambdas.pca, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
time.fimp = (proc.time() - ptm)[3]

ptm <- proc.time()
model.fimp.pg = fregression(Y ~ time | id, data, method="proximal_grad", lambda = lambdas.pca, lr = lr, thresh = 1e-4, final = "soft", maxIter = 10000, fold = 5, cv.ratio = 0.05, bins = dgrid, verbose = 0)
time.fimp.pg = (proc.time() - ptm)[3]

plot(lambdas.pca, model.fimp$loss)
mean(((ftrue - model.fimp$fit)[!is.na(simulation$fobs)])**2)
mean(((ftrue - model.fimp.pg$fit)[!is.na(simulation$fobs)])**2)

#model.X1 = fregression(X1 + X2 + Y ~ time | id, data, K=4, lambda = lambdas.pca, thresh = 1e-10, method = "fimpute", bins = dgrid, maxIter = 5000)
#model.fslr = fregression(Y ~ time + U1 + U2 + U3 + U4 | id, data, model.X1$u, K=3, thresh = 1e-10, lambda = lambdas.reg, method = "fimpute", bins = dgrid, maxIter = 5000)

# REPORT RESULTS
vv = mean((ftrue - mean(data$Y))**2)
errors = c(
  mean((ftrue - model.fpca$fit)**2) / vv,
  mean((ftrue - model.fimp$fit)**2) / vv,
  mean((ftrue - model.fimp.pg$fit)**2) / vv
  #  mean((ftrue - model.fslr$fit)**2)
)
list(errors = errors, time= c(time.fpca, time.fimp, time.fimp.pg))
}
res = mclapply(1:32, test.experiment, mc.cores = 8)

errors = c()
for (i in 1:length(res)){
  errors = rbind(errors, res[[i]]$errors)
}
times = c()
for (i in 1:length(res)){
  times = rbind(times, res[[i]]$time)
}
print(errors)
colMeans(errors)

colnames(errors) = c("fPCA",
                     "fimpute",
                     "fimpute pg")
colnames(times) = colnames(errors)

pdf("figures/simulation-errors.pdf",width=5,height=5)
boxplot(errors,ylim=c(0,1))
dev.off()
pdf("figures/simulation-times.pdf",width=5,height=5)
boxplot(times)
dev.off()

