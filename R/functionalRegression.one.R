#' Run sparse functional regression based on latent values
#'
#' @noRd
# @export
functionalRegression.one= function(Y, X, basis, lambda=0, maxIter=1e5, thresh = 1e-4, K = dim(X)[2]){
  Y = Y$train
  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = 0
  err = 1e9
  #  X = svd(X)$u
  X = cbind(1,X)

  dims = 1:K
  err = sqrt(mean( ((mean(Y,na.rm = TRUE) - Y)[!ynas])**2))

  for (i in 1:maxIter){
    Yfill[ynas] = Yhat[ynas]

    # Find B
    B = ginv(X) %*% Yfill %*% basis   # The most basic multivariate regression
    Bsvd = svd(B)
    # regularize B
    D = Bsvd$d - lambda
    D[D<0] = 0
    Dm = diag(D,length(D),length(D)) # allows D be 1x1
    B = Bsvd$u %*% Dm  %*% t(Bsvd$v)
    Yhat.new = X %*% B %*% t(basis)   #Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims]))

    # Regularize Y?
    # Ysvd = svd(Yfill %*% basis)
    # B = ginv(X) %*% Ysvd$u   # The most basic multivariate regression
    # Yu = X %*% B   #Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims]))
    # D = Ysvd$d - lambda
    # D[D<0] = 0
    # Yhat.new = Yu %*% diag(D) %*% t(Ysvd$v)  %*% t(basis)

    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-15)
    if (ratio < thresh){
      break
    }
    Yhat = Yhat.new

    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
  }
  list(fit = Yhat, coef = B, id = row.names(Y), grid = as.numeric(colnames(Y)), err = err)
}

