#' @export
functionalImpute = function(Y, basis, K = ncol(basis), maxIter = 10e3, thresh = 10e-4, lambda = 0){
  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = mean(Y, na.rm=TRUE) + rnorm(length(Y[]))

  dims = 1:K

  for (i in 1:maxIter){
    Yfill[ynas] = Yhat[ynas]
    Ysvd = svd(Yfill %*% basis)
    D = Ysvd$d[dims]- lambda
    D[D<0] = 0
    Yhat.new = Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims])) %*% t(basis)
    ratio = norm(Yhat.new - Yhat,"F") / norm(Yhat,type = "F")
    if (ratio < thresh)
     break
    Yhat = Yhat.new

    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
    cat(err,sum(D > 0),"\n")
  }
  list(fit=Yhat, D=D, U=Ysvd$u[,dims], v=(t(Ysvd$v[,dims])) %*% t(basis))
}
