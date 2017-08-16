#' Run sparse functional regression based on latent values
#'
#' @noRd
# @export
functionalRegression = function(Y, X, basis, lambda=0, maxIter=1e5, thresh = 1e-4, K = dim(X)[2]){
  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = 0
  err = 1e9
  #  X = svd(X)$u

  dims = 1:K
  err = 0
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

    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-10)
    if (ratio < thresh){
      break
    }
    Yhat = Yhat.new

    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
    cat(err,"\n")
  }
  list(fit = Yhat, id = row.names(Y), grid = as.numeric(colnames(Y)))
}

#' Parse the formula "response ~ covariates | groups"
#' to lists: response, covariates, groups
#'
#' syntax var1:var2 will add 2 variables to a corresponding list
#'
#' @noRd
# @export
parse.formula <- function(formula) {
  vars <- terms(as.formula(formula))
  y <- if(attr(vars, "response"))
    nlme::getResponseFormula(formula)
  x <- nlme::getCovariateFormula(formula)
  z <- nlme::getGroupsFormula(formula)
  list(response = all.vars(y),
       covariates = all.vars(x),
       groups = all.vars(z))
}
