functionalImpute.one = function(Y, basis, K, maxIter, thresh, lambda){
  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = 0
  err = 1e9

  dims = 1:K
  for (i in 1:maxIter){
    Yfill[ynas] = Yhat[ynas]
    Ysvd = svd(Yfill %*% basis)
    D = Ysvd$d[dims]- lambda
    D[D<0] = 0
    Yhat.new = Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims])) %*% t(basis)
    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 10e-10)
    if (ratio < thresh)
      break
    Yhat = Yhat.new

    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
    # cat(err,sum(D > 0),"\n")
  }
  list(fit=Yhat, D=D, U=Ysvd$u[,dims], v=(t(Ysvd$v[,dims])) %*% t(basis), err=err, lambda = lambda)

}

#' @export
functionalImpute = function(Y, basis = fc.basis(), K = ncol(basis), maxIter = 10e3, thresh = 10e-4, lambda = 0){
  err = 1e9
  best = NULL
  bestK = K

  smpl = fc.sample(Y)
  for (l in lambda)
  {
    model = functionalImpute.one(smpl$train, basis, K, maxIter, thresh, l)
    err.new = sqrt(mean((smp$test - model$fit)[smp$test.mask]**2))
    if (err.new < err){
      err = err.new
      best = model
      bestLambda = l
      bestK = sum(model$D > 1e-10)
    }
  }
  functionalImpute.one(Y, basis, bestK, maxIter, thresh, bestLambda)
}
