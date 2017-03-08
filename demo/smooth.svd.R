#install.packages("lbfgs")
library("lbfg")
load("mcomplete.Rda")
nobs = dim(X)[1]
ndim = dim(X)[2]
ncomp = 2

## TEST CASE
#U1 = matrix(rnorm(nobs*ncomp),nobs)
#V1 = matrix(rnorm(ncomp*ndim),ndim)
#X = U1%*%t(V1)

U = matrix(rnorm(nobs*ncomp),nobs)
V = matrix(rnorm(ncomp*ndim),ndim)

tovector = function(U,V)
{
  c(U,V)
}
fromvector = function(v)
{
  U = matrix(v[1:(nobs*ncomp)],nobs) 
  V = matrix(v[nobs*ncomp + 1:(ncomp*ndim)],ndim) 
  list(U=U,V=V)
}

lambda = 1

objective = function(v)
{
  R = fromvector(v)
  U = R$U
  V = R$V
  X.est = U %*% t(V)
  domain = !is.na(X)
  obj = mean((X[domain] - X.est[domain])**2) #+ lambda*sum((V[2:ndim,] - V[2:ndim - 1,])**2)
  print(obj)
  obj
}

gradient = function(v)
{
  R = fromvector(v)
  U = R$U
  V = R$V
  X.est = U %*% t(V)
  domain = !is.na(X)
  
  full = U%*%t(V)
  tt = X - U%*%t(V)
  full[!is.na(tt)] = tt[!is.na(tt)]
  
  Up = - full %*% V
  Vp = t(- t(U) %*% (full) )
  grad = c(Up,Vp)
#  print(grad)
  grad/length(grad)
  #sum(2*(X[domain] - X.est[domain])) 
  #+ lambda*sum((V[2:ndim,] - V[2:ndim - 1,])**2)
}

#myoptim = function(par, fn, gr){
#  for (i in 1:10000){
#    grad = gr(par)
#    par = par - grad * 0.01
#    print(fn(par))
#  }
#  list(par = par)
#}

res = optim(tovector(U,V), objective, gradient, method="L-BFGS-B")
R = fromvector(res$par)
plot(R$U[,2],U1[,2])


#######################
load("mcomplete.Rda")
X = X[,1:35]
mu = colMeans(X,na.rm = TRUE)
X = t(t(X) - colMeans(X,na.rm = TRUE))
sum.nnzero = X[1,]%*%t(X[1,]) 
sum.nnzero[] = 0
cv.est = sum.nnzero
for (r in 1:nrow(X)){
  cvv = X[r,]%*%t(X[r,]) 
  cv.est[!is.na(cvv)] = cv.est[!is.na(cvv)] + cvv[!is.na(cvv)]
  nnzero = !is.na(cvv)
  sum.nnzero = sum.nnzero + nnzero
}
Sigma = cv.est / sum.nnzero

library("MASS")
md = softImpute(Sigma)
Sigma = complete(Sigma,md)
E = eigen(Sigma)
SS = (E$vectors %*% diag(E$values)) %*% t(E$vectors)
eigen(SS)

matplot(t(mvrnorm(30, mu, SS)), t='l')
