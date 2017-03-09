# library("devtools")
# devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
# devtools::install(".")

library("fpca")
library("fcomplete")
library("fda")

d = 7
K = 2
dgrid = 501

# Set up basis
basis = fda::create.fourier.basis(c(0,1), d)
#basis = fda::create.bspline.basis(c(0,1), d, 4)
S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)

# Generate some cool positive definite matrix
V = svd(matrix(rnorm(d*d),d))$v
D = diag(c(1,0.9,0.5,exp(-(3:(d-1))))) * 500
Sigma = V %*% D %*% t(V)
Ycoef = mvrnorm(n = 1000, Sigma = Sigma, mu = rep(0,d))
Ytrue = Ycoef %*% t(S)

# Remove 90% of points
nel = prod(dim(Ytrue))
nna = ceiling(nel*0.99)
SigmaBig = genPositiveDefMat(dgrid)$Sigma
noise = mvrnorm(n = 1000, SigmaBig, mu = rep(0,dgrid)) * 0.5
Ynoise = Ytrue+noise
Y = Ynoise
par(mfrow=c(1,1))
matplot(t(Ytrue)[,1:10],t='l')
matplot(t(Y)[,1:10],t='l')
Y[sample(nel)[1:nna]] = NA

# Impute functions
Yhat = functionalImpute(Y, S, 3, maxIter = 250)
mean((Ytrue[idx,] - (Yhat[idx,]))^2) / mean(Ytrue^2)
matplot(t(Yhat)[,1:10],t='l')

# See unbelievable results

nonna = which(!is.na(Y),arr.ind = T)
nonna = nonna[order(nonna[,1]),]
vals = Y[nonna]
R = cbind(nonna, vals)
R[,2] = R[,2] / dgrid
R = R[,c(1,3,2)]

# model and predict
M = d + 2
K = 2
model = fpca.mle(R, M, K, ini.method="EM")
fpcs = fpca.score(R, model$grid, model$fitted_mean, model$eigenvalues, model$eigenfunctions, model$error_var, K)
pred = fpca.pred(fpcs, model$fitted_mean, model$eigenfunctions)

# plot results
matplot(pred, t='l')
idx = unique(R[,1])
Ycmp = Y[idx,]

tocmp = !is.na(Ycmp)
mean((Ycmp[tocmp] - t(pred)[tocmp])^2) / mean((Ycmp[tocmp])^2)
mean((Ycmp[tocmp] - (Yhat[idx,])[tocmp])^2) / mean((Ycmp[tocmp])^2)

mean((Ytrue[idx,] - t(pred))^2) / mean(Ytrue^2)
mean((Ytrue[idx,] - (Yhat[idx,]))^2) / mean(Ytrue^2)

# Plots
par(mfrow=c(2,3))
matplot(t(Ytrue[idx,])[,1:10],t='l',ylim = c(min(Ynoise),max(Ynoise)))
matplot(t(Ynoise[idx,])[,1:10],t='l',ylim = c(min(Ynoise),max(Ynoise)))
matplot(t(Y[idx,])[,1:10],t='p',pch=20,ylim = c(min(Ynoise),max(Ynoise)))
matplot(t(Yhat[idx,])[,1:10],t='l',ylim = c(min(Ynoise),max(Ynoise)))
matplot((pred)[,1:10],t='l',ylim = c(min(Ynoise),max(Ynoise)))
frame()
