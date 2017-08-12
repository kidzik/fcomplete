#' @export
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
    D = Bsvd$d - lambda
    D[D<0] = 0
    Dm = diag(D,length(D),length(D)) # allows D be 1x1
    B = Bsvd$u %*% Dm  %*% t(Bsvd$v)
    Yhat.new = X %*% B %*% t(basis)   #Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims]))

    # Regularize?
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

#' @export
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

#' @export
fregression = function(formula, data, bins = 51, method = c("fimpute", "fpca"), lambda = c(0), lambda.reg = 0, K = NULL, K.reg = NULL, thresh = 1e-5)
{
  if (length(method) > 1)
    method = "fimpute"

  d = 7
  if (is.null(K))
    K = d
  basis = fc.basis(d = d, dgrid = bins)

  vars = parse.formula(formula)
  subj.var = vars$groups

  if (length(vars$response) == 1)
  {
    # If no Y and we do unsupervised learning
    time.var = vars$response[1]
  }
  else
  {
    # If there is Y to regress
    y.var = vars$response[1]
    time.var = vars$response[2]
    Y = na.omit(data[,c(subj.var, time.var, y.var)])
    Y.wide = fc.long2wide(Y[,1], as.numeric(Y[,2]), as.numeric(Y[,3]), bins = bins)
  }

  # Case 1: Y ~ 1 -- do functional impute
  if (length(vars$covariates) == 0)
  {
    if (method == "fpcs"){
      return (fc.fpca(Y, d = d, K = K, grid.l = 0:(bins-1)/(bins-1)))
    }
    else {
      return (functionalMultiImpute(Y.wide, basis = basis, lambda = lambda, K = K, thresh = thresh))
    }
  }

  # Assome there are covariates
  X.long = list()
  X.wide = list()

  nvars = length(vars$covariates)
  for (i in 1:nvars)
  {
    X.long[[i]] = na.omit(data[,c(subj.var, time.var, vars$covariates[i])])
    X.wide[[i]] = fc.long2wide(X.long[[i]][,1], as.numeric(X.long[[i]][,2]), as.numeric(X.long[[i]][,3]), bins = bins)
  }

  # Case 2: 1 ~ X -- do unsupervised learning
  if (length(vars$response) == 1)
  {
    args = X.wide
    args$basis = basis
    args$lambda = lambda
    args$thresh = thresh
    return(do.call(functionalMultiImpute, args))
  }

  # Case 3: Y ~ X -- do regression
  models = list()
  combinedU = c()
  for (i in 1:nvars)
  {
    if (method == "fpcs"){
      models[[i]] = fc.fpca(X.long[[i]], d = d, K = 3, grid.l = 0:(bins-1)/(bins-1))
      combinedU = cbind(combinedU, models[[i]]$fpcs)
    }
    else {
      models[[i]] = functionalMultiImpute(X.wide[[i]], basis = basis, lambda = lambda, thresh = thresh, K = K)
      combinedU = cbind(combinedU, models[[i]]$u)
    }
  }

  if (is.null(K.reg))
    K.reg = ncol(Y.wide)
  functionalRegression(Y.wide, combinedU, basis, lambda = lambda.reg, K = K.reg, thresh = thresh)
}
