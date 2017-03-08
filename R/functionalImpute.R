functionalImpute = function(Y, basis, K = NULL, maxIter = 100){
  basis = svd(basis)$u
  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = mean(Y, na.rm=TRUE)
  if (is.null(K))
    K = ncol(basis)
  dims = 1:K

  for (i in 1:maxIter){
    Yfill[ynas] = Yhat[ynas]
    Ysvd = svd(Yfill %*% basis)
    Yhat = Ysvd$u[,dims] %*% (Ysvd$d[dims] * t(Ysvd$v[,dims])) %*% t(basis)
  }
  Yhat
}
