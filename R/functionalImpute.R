#' @export
functionalImpute = function(Y, basis, K = NULL, maxIter = 1000, lambda = 0){
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
    D = Ysvd$d[dims]- lambda
    D[D<0] = 0
    Yhat.new = Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims])) %*% t(basis)
    ratio = norm(Yhat.new - Yhat,"F") / norm(Yhat,type = "F")
    # if (i %% 10 == 0)
    #   print(ratio)
    # if (ratio < 1e-5)
    #   break
    Yhat = Yhat.new
  }
  list(fit=Yhat, D=D, U=Ysvd$u[,dims], v=(t(Ysvd$v[,dims])) %*% t(basis))
}
