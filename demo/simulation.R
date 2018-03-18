library("devtools")
# devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
# devtools::install(".")

library("fpca")
library("fcomplete")
library("clusterGeneration")
library("fda")
library("fMultivar")
#source("R/mixedmodel.R")

ntrials = 1
RES = array(0,c(2,2,ntrials))

set.seed(120)
d = 9
dgrid = 51
n = 150
missing = 0.05
grid = 0:(dgrid-1)/(dgrid-1)

# Set up basis
#basis = fda::create.fourier.basis(c(0,1), d)
basis = fda::create.bspline.basis(c(0,1), d, 3)
Bn = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)
B = svd(Bn)$u

# Generate some cool positive definite matrix
V = svd(matrix(rnorm(d*d),d))$v
D = diag(c(1,0.9,0.5,exp(-(3:(d-1))))) * 100
Sigma = V %*% D %*% t(V)

noise = mvrnorm(n = n, diag(d), mu = rep(0,d)) * 5
Ycoef = mvrnorm(n = n, Sigma = Sigma, mu = rnorm(d) * 0)
Ytrue = (Ycoef + noise) %*% t(B)

# Remove 90% of points
nel = prod(dim(Ytrue))
nna = floor(nel * missing)
noise = mvrnorm(n = n, diag(dgrid), mu = rep(0,dgrid)) * 0
Ynoise = Ytrue+noise
Y = Ynoise
par(mfrow=c(1,1))
# matplot(t(Ytrue)[,1:10],t='l')
matplot(t(Y)[,1:10],t='l')
# scores = mm.scores(M, Y)
# Y[sample(nel)[1:nna]] = NA

library("psych")
facta = fa(Y %*%B, 3, scores = "regression", rotate = "none", covar = TRUE)
R = factor.scores(Y %*% B, facta)
Yhat = R$scores %*% t(facta$loadings) %*% t(B)
mean((Ytrue - (Yhat))^2) / mean(Ytrue^2)

# Impute functions
fimodel = functionalImpute(Y, B, 3, lambda = 5)
Yhat = fimodel$fit
mean((Ytrue - (Yhat))^2) / mean(Ytrue^2)
# matplot(t(Yhat)[,1:10],t='l')

# See unbelievable results

nonna = which(!is.na(Y),arr.ind = T)
nonna = nonna[order(nonna[,1]),]
vals = Y[nonna]
R = cbind(nonna, vals)
R[,2] = R[,2] / dgrid
R = R[,c(1,3,2)]

# model and predict
M = d + 4
K = 4
model = fpca.mle(R, M, K, ini.method="EM",grids = grid)
fpcs = fpca.score(R, model$grid, model$fitted_mean, model$eigenvalues, model$eigenfunctions, model$error_var, K)
pred = fpca.pred(fpcs, model$fitted_mean, model$eigenfunctions)

# plot results
matplot(pred, t='l')
idx = unique(R[,1])
Ycmp = Y[idx,]

tocmp = !is.na(Ycmp)

# NMSE OBSERVED
RES[1,1,i] = mean((Ycmp[tocmp] - (Yhat[idx,])[tocmp])^2) / mean((Ycmp[tocmp])^2)
RES[2,1,i] = mean((Ycmp[tocmp] - t(pred)[tocmp])^2) / mean((Ycmp[tocmp])^2)

# NMSE TRUE
RES[1,2,i] = mean((Ytrue[idx,] - (Yhat[idx,]))^2) / mean(Ytrue^2)
RES[2,2,i] = mean((Ytrue[idx,] - t(pred))^2) / mean(Ytrue^2)

apply(RES,1:2,mean)
apply(RES,1:2,sd)

# Plots
## Figure 1: Spline basis and it's orthonormalized version
par(mfrow=c(1,2), cex=1.5)
matplot(grid,S,t='l',ylab="value")
title("Spline basis")
matplot(grid,svd(S)$u,t='l',ylab="value")
title("Orthonormalized spline basis")

## Figure 2: Data and fits
par(mfrow=c(2,2), cex=1.5)
bd = c(min(Ynoise),max(Ynoise))
bd = c(-2,2)

matplot(grid, t(Ytrue[idx,])[,1:5],t='l',ylim = bd, ylab="value")
title("TRUE")
matplot(grid, t(Ynoise[idx,])[,1:5],t='l',ylim = bd* 3, ylab="value")
title("TRUE + NOISE")
# matplot(grid, t(Y[idx,])[,1:5],t='p',pch=20,ylim = bd, ylab="value")
# title("SAMPLED TRUE + NOISE")
matplot(grid, t(Yhat[idx,])[,1:5],t='l',ylim = bd, ylab="value")
title("SOFT-IMPUTE")
matplot(grid, (pred)[,1:5],t='l',ylim = bd, ylab="value")
title("SPARSE FUNC.")

## Figure 3: Comparison of components
par(mfrow=c(1,2), cex=1.5)
lens = apply(model$eigenfunctions,1,function(x){ sqrt(sum(x^2))})

matplot(t(diag(c(1,-1)) %*% fimodel$v[1:3,]),t='l')
legend("topleft", inset=.02, legend=c("1st component", "2nd component"),
       col=c("black", "red"), lty=1:2, cex=1)
title("Soft-impute")

matplot(t(diag(1/lens) %*% model$eigenfunctions)[,1:3],t='l')
legend("topleft", inset=.02, legend=c("1st component", "2nd component"),
       col=c("black", "red"), lty=1:2, cex=1)
title("Sparse fPCA")

## Figure 4: Comparison of scores
par(mfrow=c(1,2), cex=1.5)
sparse.coeffs = fpcs %*% diag(lens)
svd.coeffs = fimodel$U %*% diag(fimodel$D * c(1,-1,1,1,1))

plot(sparse.coeffs[,1],t='l',ylab = "score")
lines(svd.coeffs[,1],t='l',col=2)
legend("topright", inset=.02, legend=c("Sparse fPCA", "Soft-impute"),
       col=c("black", "red"), lty=1, cex=1)
title("1st score")

plot(sparse.coeffs[,2],t='l',ylab = "score")
lines(svd.coeffs[,2],t='l',col=2)
legend("topright", inset=.02, legend=c("Sparse fPCA", "Soft-impute"),
       col=c("black", "red"), lty=1, cex=1)
title("2nd score")

