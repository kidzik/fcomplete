library("devtools")
#devtools::install_git("https://github.com/kidzik/fpca.git",reload = TRUE)
devtools::install(".")

library("clusterGeneration")
library("fpca")
library("fcomplete")
library("fda")
library("mvtnorm")

noise_mag = 0.2
d = 7
K = 3
dgrid = 51
clean = 0.85
n = 2000

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
dim(wide.Z)

devtools::install(".") ; library("fcomplete")

long.train = fc.wide2long(smp.Y$train)
fpca.model.Y = fc.fpca(long.train[,],d = d,K=c(K),grid.l = 0:50/50)

long.train = fc.wide2long(smp.X$train)
fpca.model.X = fc.fpca(long.train[,],d = d,K=c(K),grid.l = 0:50/50)

long.train = fc.wide2long(smp.Z$train)
fpca.model.Z = fc.fpca(long.train[,],d = d,K=c(K),grid.l = 0:50/50)

functionalRegression = function(Y, basis, U, lambda=0, maxIter=1e5, thresh = 1e-5, K = dim(U)[2]){
  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = 0
  err = 1e9
#  U = svd(U)$u

  dims = 1:K
  err = 0
  for (i in 1:maxIter){
    Yfill[ynas] = Yhat[ynas]

    # Find B
    B = ginv(U) %*% Yfill %*% basis   # The most basic multivariate regression
    Bsvd = svd(B)
    D = Bsvd$d - lambda
    D[D<0] = 0
    B = Bsvd$u %*% diag(D) %*% t(Bsvd$v)
    Yhat.new = U %*% B %*% t(basis)   #Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims]))

    # Regularize?
    # Ysvd = svd(Yfill %*% basis)
    # B = ginv(U) %*% Ysvd$u   # The most basic multivariate regression
    # Yu = U %*% B   #Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims]))
    # D = Ysvd$d - lambda
    # D[D<0] = 0
    # Yhat.new = Yu %*% diag(D) %*% t(Ysvd$v)  %*% t(basis)

    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-10)
    if (ratio < thresh){
      break
    }
    Yhat = Yhat.new

    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
    cat(err,"\n")
  }
  list(fit=Yhat)

}

dcmp = 5
U = cbind(fpca.model.X$fpcs,fpca.model.Z$fpcs)
freg = functionalRegression(Y, fc.basis(d, "splines", dgrid = dgrid), U, thresh = 1e-5, lambda=0.1)
matplot(t(smp.Y$train[ind,]),t='p',pch = 'X')
matplot(t(freg$fit[ind,]),t='l',add=T,lwd=3,lty=1)
matplot(t(fpca.model.Y$fit[ind,]),t='l',add=T,lwd=3,lty=2)
matplot(t(smp.Y$test[ind,]),t='p',pch = 'O',add=T)
sqrt(mean((Ytrue[idx,] - freg$fit)**2))

devtools::install(".") ; library("fcomplete")
# fpca.model$sigma.est * fpca.model$sigma.est * sigma.factors

scales = c(2,3.5,5,7.5,10)
func.impute = functionalMultiImpute(smp.Y$train,
                                    smp.X$train * 3,
                                    smp.Z$train * 3,
                                    basis = fc.basis(d, "splines", dgrid = dgrid),
                                    maxIter = 10e5, thresh= 1e-5,
                                    lambda = fpca.model.X$sigma.est * fpca.model.X$sigma.est * scales,
                                    K = fpca.model$selected_model[2],
                                    final="soft")
#scales = c(0.1, 0.3, 0.5, 0.7, 1, 1.3)
#func.impute = functionalMultiImpute(smp.Y$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-6, lambda = fpca.model$sigma.est * fpca.model$sigma.est * scales,K=3,final="soft")
#func.impute = functionalMultiImpute(smp.Y$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-10, lambda = fpca.model$sigma.est * fpca.model$sigma.est * c(1,1.5,2,3),K=4)
#func.impute = functionalMultiImpute(smp.Y$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-4, lambda = fpca.model$sigma.est * fpca.model$sigma.est * c(0.1,0.3,0.7,1,1.5,2,3), K=2)
#func.impute = functionalMultiImpute.one(smp.Y, smp.X, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-5, lambda = fpca.model$sigma.est * fpca.model$sigma.est * c(1.2), K=2)
func.impute$fit = func.impute$fit[[1]]

#func.impute = functionalImpute.one(smp.Y$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-5, lambda = fpca.model$sigma.est * fpca.model$sigma.est *1.5, K=2)
#func.impute$fit = func.impute$fit[[1]]

#func.impute = functionalMultiImpute(smp.Y$train, basis = fc.basis(d, "splines", dgrid = dgrid), maxIter = 10e5, thresh= 1e-3, lambda = fpca.model$sigma.est * fpca.model$sigma.est * sigma.factors)
mean.impute = fc.mean(smp.Y$train)

par(cex=2)
ind = 1:3 + 20
matplot(t(smp.Y$train[ind,]),t='p',pch = 'X')
matplot(t(mean.impute[ind,]),t='l',add=T,lwd=3)
matplot(t(smp.Y$test[ind,]),t='p',pch = 'O',add=T)

matplot(t(smp.Y$train[ind,]),t='p',pch = 'X')
matplot(t(func.impute$fit[ind,]),t='l',add=T,lwd=3,lty=2)
matplot(t(smp.Y$test[ind,]),t='p',pch = 'O',add=T)

matplot(t(smp.Y$train[ind,]),t='p',pch = 'X')
matplot(t(fpca.model.X$fit[ind,]),t='l',add=T,lwd=3,lty=2)
matplot(t(smp.Y$test[ind,]),t='p',pch = 'O',add=T)

ensamble = (fpca.model.X$fit + func.impute$fit)/2
matplot(t(ensamble[ind,]),t='l',lwd=3)
matplot(t(smp.Y$train[ind,]),t='p',pch = 'X',add=T)
matplot(t(smp.Y$test[ind,]),t='p',pch = 'O',add=T)

m1 = sqrt(mean((smp.Y$test - func.impute$fit)[smp.Y$test.mask]**2))
m0 = sqrt(mean((smp.Y$test - mean.impute)[smp.Y$test.mask]**2))
m2 = sqrt(mean((smp.Y$test - fpca.model.X$fit)[smp.Y$test.mask]**2))
m3 = sqrt(mean((smp.Y$test - ensamble )[smp.Y$test.mask]**2))
m4 = sqrt(mean((smp.Y$test - freg$fit)[smp.Y$test.mask]**2))

tbl.smp = rbind(
  c(m0,100*(1-m0/m0)),
  c(m1,100*(1-m1/m0)),
  c(m2,100*(1-m2/m0)),
  c(m3,100*(1-m3/m0)),
  c(m4,100*(1-m4/m0))
)
colnames(tbl.smp) = c("MSE","% expl")

#cat("SAMPLE:\nmean impute:\t",m0,"\t",0,"%\nours:\t\t",m1,"\t",(100*m1/m0),"%\nfpca:\t\t",m2,"\t",(100*m2/m0),"%\nensamble:\t",m3,"\t",(100*m3/m0),"%")

m1 = sqrt(mean((Ytrue[idx,] - func.impute$fit)**2))
m0 = sqrt(mean((Ytrue[idx,] - mean.impute)**2))
m2 = sqrt(mean((Ytrue[idx,] - fpca.model.X$fit)**2))
m3 = sqrt(mean((Ytrue[idx,] - ensamble )**2))
m4 = sqrt(mean((Ytrue[idx,] - freg$fit)**2))

tbl.true = rbind(
  c(m0,100*(1-m0/m0)),
  c(m1,100*(1-m1/m0)),
  c(m2,100*(1-m2/m0)),
  c(m3,100*(1-m3/m0)),
  c(m4,100*(1-m4/m0))
)
colnames(tbl.true) = c("MSE","% expl")

tbl.final = cbind(c(0, fpca.model.X$selected_model[2], sum(func.impute$D > 0.01), NA, NA ), tbl.smp, tbl.true)
colnames(tbl.final)[1] = "K"
rownames(tbl.final) = c("mean","ours","fpca","ensamble","regression")

#cat("TRUE:\nmean impute:\t",m0,"\t",0,"%\nours:\t\t",m1,"\t",(100*m1/m0),"%\nfpca:\t\t",m2,"\t",(100*m2/m0),"%\nensamble:\t",m3,"\t",(100*m3/m0),"%")
print(tbl.final)

matplot(t(Ytrue[idx,][ind,]),t='l',lty=1,lwd=4)
matplot(t(func.impute$fit[ind,]),t='l',lty=2,add=T,lwd=2)

matplot(t(Ytrue[idx,][ind,]),t='l',lty=1,lwd=4)
matplot(t(fpca.model$fit[ind,]),t='l',lty=2, add=T,lwd=2)

matplot(t(Ytrue[idx,][ind,]),t='l',lty=1,lwd=4)
matplot(t(freg$fit[ind,]),t='l',lty=2, add=T,lwd=2)

