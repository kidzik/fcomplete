# library("devtools")
# devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
# devtools::install(".")

library("clusterGeneration")
library("fpca")
library("fcomplete")
library("fda")

d = 7
K = 2
noise_mag = 0.5
dgrid = 501
clean = 0.995

# Set up a basis
basis = fda::create.fourier.basis(c(0,1), d)
#basis = fda::create.bspline.basis(c(0,1), d, 4)
S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)

# Generate some cool positive definite matrix
V = svd(matrix(rnorm(d*d),d))$v
D = diag(c(1,0.9,0.5,exp(-(3:(d-1))))) * 500
Sigma = V %*% D %*% t(V)
Ycoef = mvrnorm(n = 1000, Sigma = Sigma, mu = rep(0,d))
Ytrue = Ycoef %*% t(S)

# Remove 99% of points
nel = prod(dim(Ytrue))
nna = ceiling(nel*clean)
SigmaBig = genPositiveDefMat(dgrid)$Sigma
noise = mvrnorm(n = 1000, SigmaBig, mu = rep(0,dgrid)) * noise_mag
Ynoise = Ytrue+noise
Y = Ynoise
Y[sample(nel)[1:nna]] = NA

# model and predict
Ms = 5:10
Ks = 2:4
errs.fpca = matrix(NA,length(Ks),length(Ms))
colnames(errs.fpca) = Ms
rownames(errs.fpca) = Ks

nonna = which(!is.na(Y),arr.ind = T)
nonna = nonna[order(nonna[,1]),]
vals = Y[nonna]
R = cbind(nonna, vals)
R[,2] = R[,2] / dgrid
R = R[,c(1,3,2)]
idx = unique(R[,1])

for (i in 1:length(Ks)){
  for (j in 1:length(Ms)){
    K = Ks[i]
    M = Ms[j]
    model = fpca.mle(R, M, K, ini.method="EM")
    fpcs = fpca.score(R, model$grid, model$fitted_mean, model$eigenvalues, model$eigenfunctions, model$error_var, K)
    pred = fpca.pred(fpcs, model$fitted_mean, model$eigenfunctions)
    errs.fpca[i,j] = mean((Ytrue[idx,] - t(pred))^2) / mean(Ytrue^2)
    print(errs.fpca[i,j])
  }
}

# Impute functions
Ks = 1:7
iters = 1:4 * 50
errs = matrix(NA,length(Ks),length(iters))
colnames(errs) = iters
rownames(errs) = Ks
for (i in 1:length(iters)){
  print(iters[i])
  maxIter = iters[i]
  for (K in Ks){
    Yhat = functionalImpute(Y, S, K, maxIter = maxIter)
    errs[K,i] = mean((Ytrue[idx,] - (Yhat[idx,]))^2) / mean(Ytrue^2)
    print(errs[K,i])
  }
}

