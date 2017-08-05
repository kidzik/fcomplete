library("devtools")
#devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
devtools::install(".")

library("clusterGeneration")
library("fpca")
library("fcomplete")
library("fda")
library("mvtnorm")

noise_mag = 0.1
d = 5
K = 3
dgrid = 51
clean = 0.8
n = 1000

# Set up a basis
#basis = fda::create.fourier.basis(c(0,1), d)
basis = fda::create.bspline.basis(c(0,1), d, 4)
S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)

# Generate some cool positive definite matrix
generate.matrix = function(){
  V = svd(matrix(rnorm(d*d),d))$v
  D = diag(c(1,0.9,0.5,exp(-(3:(d-1))))) * 500
  Sigma1 = V %*% D %*% t(V)
  V = svd(matrix(rnorm(d*d),d))$v
  D = diag(c(1.3,0.4,0.4,exp(-(3:(d-1))))) * 500
  Sigma2 = V %*% D %*% t(V)
  Ycoef1 = rmvnorm(n = n, sigma = Sigma1, mean = rnorm(d))
  Ycoef2 = rmvnorm(n = n, sigma = Sigma2, mean = rnorm(d))
  subst = runif(n) > 0.2
  Ycoef = Ycoef1
  Ycoef[subst,] = Ycoef2[subst,]
  Ycoef
}
B = matrix(rnorm(d**2),d,d)
Xcoef = generate.matrix()
Zcoef = generate.matrix()
Ycoef = Xcoef + Zcoef #+ generate.matrix() * 0.1

Ztrue = Zcoef %*% t(S)
Ytrue = Ycoef %*% t(S)
Xtrue = Xcoef %*% t(S)

# Remove 99% of points
nel = prod(dim(Ytrue))
nna = ceiling(nel*clean)
SigmaBig = genPositiveDefMat(dgrid)$Sigma
noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag * 10
Ynoise = Ytrue+noise
noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag
Xnoise = Xtrue+noise
noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag
Znoise = Ztrue+noise

Y = Ynoise
remove.points = sample(nel)[1:nna]
Y[remove.points] = NA
X = Xnoise
X[remove.points] = NA
Z = Znoise
Z[remove.points] = NA

# model and predict
Ms = d #5:10
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

wide.Y = fc.long2wide(R[,1],R[,3],R[,2],bins = dgrid)
smp.Y = fc.sample(wide.Y)

nonna = which(!is.na(X),arr.ind = T)
nonna = nonna[order(nonna[,1]),]
vals = X[nonna]
R = cbind(nonna, vals)
R[,2] = R[,2] / dgrid
R = R[,c(1,3,2)]
idx = unique(R[,1])

wide.X = fc.long2wide(R[,1],R[,3],R[,2],bins = dgrid)
smp.X = fcomplete::apply.mask(wide.X, smp.Y)

nonna = which(!is.na(Z),arr.ind = T)
nonna = nonna[order(nonna[,1]),]
vals = Z[nonna]
R = cbind(nonna, vals)
R[,2] = R[,2] / dgrid
R = R[,c(1,3,2)]
idx = unique(R[,1])

wide.Z = fc.long2wide(R[,1],R[,3],R[,2],bins = dgrid)
smp.Z = fcomplete::apply.mask(wide.Z, smp.Y)

long.train = fc.wide2long(smp.Y$train)
fpca.model = fc.fpca(long.train[,],d = 7,K=c(K-1,K,K+1),grid.l = 0:50/50)

sigma.factors = c(0.1, 0.5, 0.7, 1, 1.5)

devtools::install(".")
library("fcomplete")
# fpca.model$sigma.est * fpca.model$sigma.est * sigma.factors
func.impute = functionalMultiImpute(smp.Y$train, smp.X$train, smp.Z$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-4, lambda = fpca.model$sigma.est * fpca.model$sigma.est * c(0.1,0.3,0.7,1,1.5,2,3), K=2)
#func.impute = functionalMultiImpute.one(smp.Y, smp.X, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-5, lambda = fpca.model$sigma.est * fpca.model$sigma.est * c(1.2), K=2)
func.impute$fit = func.impute$fit[[1]]

#func.impute = functionalImpute.one(smp.Y$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-5, lambda = fpca.model$sigma.est * fpca.model$sigma.est *1.5, K=2)
#func.impute$fit = func.impute$fit[[1]]

#func.impute = functionalMultiImpute(smp.Y$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-3, lambda = fpca.model$sigma.est * fpca.model$sigma.est * sigma.factors)
mean.impute = fc.mean(smp.Y$train)

ind = 10:15
matplot(t(smp.Y$train[smp.Y$test.rows[ind],]),t='p',pch = 'X')
matplot(t(mean.impute[smp.Y$test.rows[ind],]),t='l',add=T)
matplot(t(smp.Y$test[smp.Y$test.rows[ind],]),t='p',pch = 'O',add=T)

matplot(t(smp.Y$train[smp.Y$test.rows[ind],]),t='p',pch = 'X')
matplot(t(func.impute$fit[smp.Y$test.rows[ind],]),t='l',add=T)
matplot(t(smp.Y$test[smp.Y$test.rows[ind],]),t='p',pch = 'O',add=T)

matplot(t(smp.Y$train[smp.Y$test.rows[ind],]),t='p',pch = 'X')
matplot(t(fpca.model$fit[smp.Y$test.rows[ind],]),t='l',add=T)
matplot(t(smp.Y$test[smp.Y$test.rows[ind],]),t='p',pch = 'O',add=T)

ensamble = (fpca.model$fit + func.impute$fit)/2
matplot(t(ensamble[smp.Y$test.rows[ind],]),t='l')
matplot(t(smp.Y$train[smp.Y$test.rows[ind],]),t='p',pch = 'X',add=T)
matplot(t(smp.Y$test[smp.Y$test.rows[ind],]),t='p',pch = 'O',add=T)

m1 = sqrt(mean((smp.Y$test - func.impute$fit)[smp.Y$test.mask]**2))
m0 = sqrt(mean((smp.Y$test - mean.impute)[smp.Y$test.mask]**2))
m2 = sqrt(mean((smp.Y$test - fpca.model$fit)[smp.Y$test.mask]**2))
m3 = sqrt(mean((smp.Y$test - ensamble )[smp.Y$test.mask]**2))

cat("SAMPLE:\nmean impute:\t",m0,"\nours:\t\t",m1,"\nfpca:\t\t",m2,"\nensamble:\t",m3)

m1 = sqrt(mean((Ytrue[idx,] - func.impute$fit)**2))
m0 = sqrt(mean((Ytrue[idx,] - mean.impute)**2))
m2 = sqrt(mean((Ytrue[idx,] - fpca.model$fit)**2))
m3 = sqrt(mean((Ytrue[idx,] - ensamble )**2))

cat("TRUE:\nmean impute:\t",m0/m0,"\nours:\t\t",m1/m0,"\nfpca:\t\t",m2/m0,"\nensamble:\t",m3/m0)

ind = 1:3 + 20
matplot(t(Ytrue[idx,][ind,]),t='l',lty=1,lwd=2)
matplot(t(func.impute$fit[smp.Y$test.rows[ind],]),t='l',lty=2,add=T)

matplot(t(Ytrue[idx,][ind,]),t='l',lty=1,lwd=2)
matplot(t(fpca.model$fit[smp.Y$test.rows[ind],]),t='l',lty=3, add=T)
