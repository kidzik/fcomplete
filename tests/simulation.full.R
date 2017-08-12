library("devtools")
library("roxygen2")
#devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
roxygenize()
devtools::install(".") ; library("fcomplete")

library("clusterGeneration")
library("fpca")
library("fda")
library("mvtnorm")

# Simulation parameters
noise_mag = 0.2
d = 7
K = 3
dgrid = 51
clear = 0.85
n = 2000

# # Set up a basis
# basis = fc.basis(d = d, "splines", norder = 4, dgrid = dgrid)
#
# X1obj = fc.generate.one(n, basis, noise_mag)
# X2obj = fc.generate.one(n, basis, noise_mag)
#
# # Generate coefficients
# coef = X1obj$coef + X2obj$coef
# Yobj = fc.generate.one(coef, basis, noise_mag)
#
# # Remove clear*100% of coefficients
# smplX1 = fc.sample(X1obj$obs, perc = clear, one.per.row = FALSE)
# X1.wide = smplX1$train
# smplX2 = fcomplete::apply.mask(X2obj$obs, smplX1)
# X2.wide = smplX2$train
# smplY = fcomplete::apply.mask(Yobj$obs, smplX1)
# Y.wide = smplY$train
#
# subjid = as.character(1:n)
# time = as.character(0:(dgrid-1)/(dgrid-1))
#
# X1.long = fc.wide2long(X1.wide, time, subjid)
# X2.long = fc.wide2long(X2.wide, time, subjid)
# data = fc.wide2long(Y.wide, time, subjid, value="Y")
# data$X1 = X1.long$value
# data$X2 = X2.long$value
###############################################################################
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
  subst = runif(n) > 0
  Ycoef = Ycoef1
  Ycoef[subst,] = Ycoef2[subst,]
  Ycoef
}
B = matrix(rnorm(d**2),d,d)
Xcoef = generate.matrix()
Zcoef = generate.matrix()
Ycoef = Xcoef[,] + Zcoef[,] #+ generate.matrix() * 0.1

Ztrue = Zcoef %*% t(S)
Xtrue = Xcoef %*% t(S)
# Ytrue = Ycoef %*% t(S)
Ytrue = Xtrue + Ztrue

# Remove 99% of points
nel = prod(dim(Ytrue))
nna = ceiling(nel*clean)
SigmaBig = genPositiveDefMat(dgrid)$Sigma

noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag
Xnoise = Xtrue + noise
noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag
Znoise = Ztrue + noise
noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag
Ynoise = Xnoise + Znoise + noise * 4
#Ynoise = Xnoise + noise * 5

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
length(idx)
data = data.frame(id=R[,1],time = R[,3],Y = R[,2])

wide.Y = fc.long2wide(R[,1],R[,3],R[,2],bins = dgrid)
smp.Y = fc.sample(wide.Y)

nonna = which(!is.na(X),arr.ind = T)
nonna = nonna[order(nonna[,1]),]
vals = X[nonna]
R = cbind(nonna, vals)
R[,2] = R[,2] / dgrid
R = R[,c(1,3,2)]
idx = unique(R[,1])
data$X1 = R[,2]

wide.X = fc.long2wide(R[,1],R[,3],R[,2],bins = dgrid)
smp.X = fcomplete::apply.mask(wide.X, smp.Y)

nonna = which(!is.na(Z),arr.ind = T)
nonna = nonna[order(nonna[,1]),]
vals = Z[nonna]
R = cbind(nonna, vals)
R[,2] = R[,2] / dgrid
R = R[,c(1,3,2)]
idx = unique(R[,1])
data$X2 = R[,2]

wide.Z = fc.long2wide(R[,1],R[,3],R[,2],bins = dgrid)
smp.Z = fcomplete::apply.mask(wide.Z, smp.Y)
dim(wide.Z)

Yobj = list(ftrue = Ytrue)
Y.wide = wide.Y
#############################################




#par(mfrow=c(1,2))
#plot(data$X1, data$Y)
#plot(data$X2, data$Y)
#summary(lm(Y ~ X1 + X2, data))

# REGRESSION
#model.fpca = fc.fpca(data[,1:3], d = d, K = K, grid.l = 0:(dgrid-1)/(dgrid-1))
model.mean = fc.mean(Y.wide)
lambdas = model.fpca$sigma.est**2 * c(0.8,1)
model.fcompress = fregression(1:time ~ Y + X1 + X2 | id, data, lambda = lambdas, K = K)
model.freg = fregression(Y:time ~ X1 + X2 | id, data, lambda = lambdas, lambda.reg = 1, method = "fpcs", K = K)

devtools::install(".") ; library("fcomplete")
model.fimpute = fregression(Y:time ~ 1 | id, data, lambda = lambdas, K = 7, thresh = 1e-7)
model.fpca = fregression(Y:time ~ 1 | id, data, lambda = lambdas, K = c(3,4,5), thresh = 1e-7, method = "fpcs")

errors = c(sqrt(mean((Yobj$ftrue[idx,] - model.mean[idx,])**2)),
  sqrt(mean((Yobj$ftrue[idx,] - model.fpca$fit)**2)),
  sqrt(mean((Yobj$ftrue[idx,] - model.fimpute$fit[[1]])**2)),
  sqrt(mean((Yobj$ftrue[idx,] - model.fcompress$fit[[1]])**2)),
  sqrt(mean((Yobj$ftrue[idx,] - model.freg$fit)**2)))

tbl.true = cbind(
  errors,
  100*(1-errors/errors[1])
)
colnames(tbl.true) = c("MSE","% expl")
rownames(tbl.true) = c("mean","fpca","fimpute","fcompress","regression")
print(tbl.true)

# Plot some examples
par(mfrow=c(2,2))
ind = 100 + 1:3
idx = as.numeric(model.freg$id)
matplot(t(Yobj$ftrue[idx,][ind,]),t='l',lty=1,lwd=4)
matplot(t(model.mean[ind,]),t='l',lty=2,add=T,lwd=2)

matplot(t(Yobj$ftrue[idx,][ind,]),t='l',lty=1,lwd=4)
matplot(t(model.fimpute$fit[[1]][ind,]),t='l',lty=2,add=T,lwd=2)

matplot(t(Yobj$ftrue[idx,][ind,]),t='l',lty=1,lwd=4)
matplot(t(model.fpca$fit[ind,]),t='l',lty=2, add=T,lwd=2)

matplot(t(Yobj$ftrue[idx,][ind,]),t='l',lty=1,lwd=4)
matplot(t(model.freg$fit[ind,]),t='l',lty=2, add=T,lwd=2)

