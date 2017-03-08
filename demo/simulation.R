library("devtools")
# devtools::install(".",reload = TRUE)
# library("fcomplete")
library("fda")
library("clusterGeneration")

d = 21
dgrid = 100

# Set up basis
basis = fda::create.fourier.basis(c(0,1), d)
#basis = fda::create.bspline.basis(c(0,1), d, 4)
S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)

# Generate data
Sigma = genPositiveDefMat(d)$Sigma
Ycoef = mvrnorm(n = 1000, Sigma, mu = (d:1 / 3) + rnorm(d))
Ytrue = Ycoef %*% t(S)

# Remove 90% of points
nel = prod(dim(Ytrue))
nna = ceiling(nel*0.95)
SigmaBig = genPositiveDefMat(dgrid)$Sigma
noise = mvrnorm(n = 1000, SigmaBig, mu = rep(0,dgrid)) * 0.01
Ynoise = Ytrue+noise
Y = Ynoise
Y[sample(nel)[1:nna]] = NA

# Impute functions
Yhat = functionalImpute(Y, S, 3)

# See unbelievable results
par(mfrow=c(2,2))
matplot(t(Ytrue)[,7:10],t='l')
matplot(t(Ynoise)[,7:10],t='l')
matplot(t(Y)[,7:10],t='p',pch=20)
matplot(t(Yhat)[,7:10],t='l')
