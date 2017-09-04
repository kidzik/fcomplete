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

    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-15)
    if (ratio < thresh){
      break
    }
    Yhat = Yhat.new

    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
  }
  list(fit = Yhat, id = row.names(Y), grid = as.numeric(colnames(Y)), err = err)
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

#' Run sparse functional regression based on latent values
#'
#' @noRd
# @export
functionalRegression = function(Y, X, basis, lambda=0, maxIter=1e5, thresh = 1e-4, K = dim(X)[2], mask = NULL){
  err = 1e9
  bestLambda = NULL
  bestModel = NULL
  meta = NULL

  cv.err = c()
  fit.err = c()

  args.smpl = list()
  if (!is.null(mask)){
    args.smpl[["Y"]] = apply.mask(Y,mask)
  }
  else {
    args.smpl[["Y"]] = fc.sample(Y)
  }
  args.smpl[["X"]] = X
  nargs = length(args.smpl)
  args.smpl[["basis"]] = basis
  args.smpl[["K"]] = K
  args.smpl[["maxIter"]] = maxIter
  args.smpl[["thresh"]] = thresh

  if (length(lambda) > 1){
    for (l in lambda)
    {
      args.smpl[["lambda"]] = l
      model = do.call(functionalRegression.one, args.smpl)

      # print(args.smpl[["Y"]]$test[args.smpl[["Y"]]$test.mask]**2)
      # print(args.smpl[["Y"]]$train[args.smpl[["Y"]]$test.mask]**2)
      # print(model$fit[args.smpl[["Y"]]$test.mask]**2)
      # print((args.smpl[["Y"]]$test - model$fit)[args.smpl[["Y"]]$test.mask]**2 )
      err.new = sqrt(mean((args.smpl[["Y"]]$test - model$fit)[args.smpl[["Y"]]$test.mask]**2))

      cat(paste("Error with lambda=",l,"\t",err.new,"\n"))
      cv.err = c(cv.err, err.new)
      fit.err = c(fit.err, model$err)

      if (err.new < err){
        err = err.new
        bestLambda = l
        bestModel = model
      }
    }
    meta = data.frame(lambda = lambda, cv.err = cv.err, fit.err = fit.err)
  }
  else {
    bestLambda = lambda
  }

  args.smpl[["lambda"]] = bestLambda
  args.smpl[["Y"]]$train = Y
  res = do.call(functionalRegression.one, args.smpl)
  res$meta = meta
  res$err.cv = err
  res$lambda = bestLambda
  res

}
