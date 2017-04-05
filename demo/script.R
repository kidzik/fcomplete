par(mfrow=c(1,1))
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

D1 = functionalImpute(D.train,basis,K=3,maxIter=20)

matplot(t(D1[1:10,]),t='l')
print(sqrt(mean(get_errors(D1)^2)))

K = 3
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
R = R[order(R[,3]), ]
R = R[order(R[,1]), ]
model = fpca.mle(R, M, K, ini.method="EM")
fpcs = fpca.score(R, model$grid,model$fitted_mean,model$eigenvalues,model$eigenfunctions,model$error_var,K)
pred.tmp = t(fpca.pred(fpcs, model$fitted_mean, model$eigenfunctions))
#pred = t(t(pred.tmp) - colMeans(pred.tmp))
pred = pred.tmp

plotids = unique(R[,1])

plot.points = function(){
  for (i in 1:5){
    pid = plotids[i]
    X = R[R[,1] == pid,3]
    Y = R[R[,1] == pid,2]
    # if (i == 1)
    #   plot(X,Y,xlim=c(0,1), ylim = c(5,35),col=i, t='o')
    # else
    lines(X,Y,xlim=c(0,1), col=i, t='o')
  }
}

matplot(1:501 / 501,t(pred[1:5,]),t='l',ylim=c(-15,35))
title("Sparse Functional PCA")
plot.points()
matplot(1:51 / 51,t(D1[1:5,]),t='l',ylim=c(-15,35))
title("Soft Impute")
plot.points()
# matplot(t(D2[1:100,]),t='l')
# points(data.train$age*51, data.train$pc1)

errors = get_errors(D1)
sqrt(mean(errors^2))
C1 = c(C1,mean(errors^2))
errors = get_errors((pred))
sqrt(mean(errors^2))
C2 = c(C2,mean(errors^2))
