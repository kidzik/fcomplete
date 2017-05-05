# library("devtools")
# devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
# devtools::install(".")

library("fpca")
library("fcomplete")
library("clusterGeneration")
library("fda")
library("fMultivar")
source("R/mixedmodel.R")
library("zoo")

ntrials = 1
RES = array(0,c(2,2,ntrials))

set.seed(120)
d = 9
dgrid = 51
n = 1000
missing = 0.99
grid = 0:(dgrid-1)/(dgrid-1)

# Set up basis
#basis = fda::create.fourier.basis(c(0,1), d)
basis = fda::create.bspline.basis(c(0,1), d, 3)
Bn = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)
B = svd(Bn)$u

# Generate some cool positive definite matrix
V = svd(matrix(rnorm(d*d),d))$v
D = diag(c(1,0.9,0.5,0*exp(-(3:(d-1))))) * 100
Sigma = V %*% D %*% t(V)

noise = mvrnorm(n = n, diag(d), mu = rep(0,d)) * 1
Ycoef = mvrnorm(n = n, Sigma = Sigma, mu = rnorm(d) * 0)
Ytrue = (Ycoef + noise) %*% t(B)

nel = prod(dim(Ytrue))
nna = floor(nel * missing)
noise = mvrnorm(n = n, diag(dgrid), mu = rep(0,dgrid)) * 0
Ynoise = Ytrue+noise
Y = Ynoise
M = mm.fit(Y,B)

X = Y %*% B

library("psych")
facta = fa(X, 3, scores = "regression", rotate = "none", covar = TRUE)
R = factor.scores(X, facta)
R$scores %*% R$weights * t(facta$loadings)
mean(sqrt(diag(facta$r - facta$model))) # sigma estimate

# "none", "varimax", "quatimax", "promax", "oblimin", "simplimax", or "cluster"
# matplot(Y)
#
# matplot(t(Y),t='l')
# matplot(,t='l')
