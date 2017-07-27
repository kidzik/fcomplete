# library("devtools")
# devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
devtools::install(".")

library("clusterGeneration")
library("fpca")
library("fcomplete")
library("fda")
library("mvtnorm")

d = 7
K = 3
noise_mag = 1
dgrid = 501
clean = 0.995
n = 200

# Set up a basis
basis = fda::create.fourier.basis(c(0,1), d)
#basis = fda::create.bspline.basis(c(0,1), d, 4)
S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)

# Generate some cool positive definite matrix
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
#Ycoef[subst,] = Ycoef2[subst,]

Ytrue = Ycoef %*% t(S)

# Remove 99% of points
nel = prod(dim(Ytrue))
nna = ceiling(nel*clean)
SigmaBig = genPositiveDefMat(dgrid)$Sigma
noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag
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

wide = fc.long2wide(R[,1],R[,3],R[,2],bins = 100)
smp = fc.sample(wide)

long.train = fc.wide2long(smp$train)
fpca.model = fc.fpca(long.train[,],d = 7,K=c(K-1,K,K+1))

sigma.factors = c(0.1, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10)
func.impute = functionalImpute(smp$train, basis = fc.basis(d, "splines"), maxIter = 10e5, thresh= 1e-5, lambda = fpca.model$sigma.est * fpca.model$sigma.est * sigma.factors)
mean.impute = fc.mean(smp$train)

ind = 10:15
matplot(t(smp$train[smp$test.rows[ind],]),t='p',pch = 'X')
matplot(t(mean.impute[smp$test.rows[ind],]),t='l',add=T)
matplot(t(smp$test[smp$test.rows[ind],]),t='p',pch = 'O',add=T)

matplot(t(smp$train[smp$test.rows[ind],]),t='p',pch = 'X')
matplot(t(func.impute$fit[smp$test.rows[ind],]),t='l',add=T)
matplot(t(smp$test[smp$test.rows[ind],]),t='p',pch = 'O',add=T)

matplot(t(smp$train[smp$test.rows[ind],]),t='p',pch = 'X')
matplot(t(fpca.model$fit[smp$test.rows[ind],]),t='l',add=T)
matplot(t(smp$test[smp$test.rows[ind],]),t='p',pch = 'O',add=T)

ensamble = (fpca.model$fit + func.impute$fit)/2
matplot(t(ensamble[smp$test.rows[ind],]),t='l')
matplot(t(smp$train[smp$test.rows[ind],]),t='p',pch = 'X',add=T)
matplot(t(smp$test[smp$test.rows[ind],]),t='p',pch = 'O',add=T)

m1 = sqrt(mean((smp$test - func.impute$fit)[smp$test.mask]**2))
m0 = sqrt(mean((smp$test - mean.impute)[smp$test.mask]**2))
m2 = sqrt(mean((smp$test - fpca.model$fit)[smp$test.mask]**2))
m3 = sqrt(mean((smp$test - ensamble )[smp$test.mask]**2))

cat("mean impute:\t",m0,"\nours:\t\t",m1,"\nfpca:\t\t",m2,"\nensamble:\t",m3)
