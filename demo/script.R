# load("~/Dropbox/sparselong/allprojections.Rda")

get_errors = function(pred){
  errors = c()
  for (i in 1:nrow(data.test)){
    id = data.test$id[i]
    nb = ceiling(data.test$age[i] * nbins)
    val = data.test$pc1[i]

    row = which(idx == id)
    err = pred[row,nb] - val
    errors = c(errors,err)
  }
  errors
}

## PREPARE DATA
idx = unique(gaitIndex$id)
nobs = length(idx)

mmin = min(gaitIndex$age)
mmax = max(gaitIndex$age)
mapage = function(x){
  (x - mmin) / (mmax - mmin)
}

nbins = 51

data = gaitIndex
data$age = mapage(data$age)

# Get some training / test set
dlast = c(data$id[2:(nobs)] - data$id[1:(nobs-1)] > 0,FALSE)

test.idx = sample(which(dlast),70)

data.test = data[test.idx,]
data.train = data[-test.idx,]

D.train.tmp = matrix(nrow=length(idx), ncol=nbins)
for (i in 1:length(idx)){
  id = idx[i]
  ss = data.train[data.train$id == id,]
  for (v in 1:nrow(ss)){
    D.train.tmp[i,ceiling(ss$age[v] * nbins)] = ss$pc1[v]
  }
}
D.train = D.train.tmp #scale(D.train.tmp,scale = TRUE)

# data.train - sparse format train
# data.test  - sparse format test
# D.train    - matrix format train
# No need for matrix Test

## Matrix completion model

d = 6
#basis.fd = fda::create.fourier.basis(c(0,1), d)
basis.fd = fda::create.bspline.basis(c(0,1), d)
basis = fda::eval.basis(evalarg = 0:(nbins-1)/(nbins-1), basisobj = basis.fd) / sqrt(nbins)

D1 = functionalImpute(D.train,basis,K=2,maxIter=i*1)

matplot(t(D.train)[,1:100],t='p',pch=20)
matplot(t(D1[1:10,]),t='l')
print(sqrt(mean(get_errors(D1)^2)))

K = 2
M = 5

## Functiona PCA model
library("fpca")

# Fix the grid in fpca.score from the fpca library
fpca.score = function (data.m, grids.u, muhat, eigenvals, eigenfuncs, sig2hat, K)
{
  temp <- table(data.m[, 1])
  n <- length(temp)
  m.l <- as.vector(temp)
  result <- matrix(0, n, K)
  N <- length(grids.u)
  evalmat <- diag(eigenvals[1:K])
  current <- 0
  eigenfuncs.u <- t(eigenfuncs)
  data.u <- matrix(as.numeric(as.vector(data.m[, -1])), nrow = nrow(data.m[,-1]), ncol = ncol(data.m[, -1]))
  for (i in 1:n) {
    Y <- as.vector(data.u[(current + 1):(current + m.l[i]), 1])
    meastime <- data.u[(current + 1):(current + m.l[i]), 2]
    gridtime <- ceiling(N * meastime)
    gridtime[gridtime <= 0] = 1 # QUICKFIX
    muy <- muhat[gridtime]
    Phiy <- matrix(eigenfuncs.u[gridtime, 1:K], ncol = K)
    Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
    temp.y <- matrix(Y - muy)
    result[i, ] <- evalmat %*% t(Phiy) %*% solve(Sigy, temp.y)
    current <- current + m.l[i]
  }
  return(result)
}

R = as.matrix(data.train)
model = fpca.mle(R, M, K, ini.method="EM")
fpcs = fpca.score(R, model$grid,model$fitted_mean,model$eigenvalues,model$eigenfunctions,model$error_var,K)
pred.tmp = t(fpca.pred(fpcs, model$fitted_mean, model$eigenfunctions))
#pred = t(t(pred.tmp) - colMeans(pred.tmp))
pred = pred.tmp

matplot(t(pred[1:10,]),t='l')
# matplot(t(D2[1:100,]),t='l')
# points(data.train$age*51, data.train$pc1)

errors = get_errors(D1)
sqrt(mean(errors^2))
errors = get_errors((pred))
sqrt(mean(errors^2))
