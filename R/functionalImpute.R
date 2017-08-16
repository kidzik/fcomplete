#' Project \code{Y} potentially composed of multiple variables \code{X1,X2,...,Xp}
#' onto \code{basis^p} space
#' @noRd
project.on.basis = function(Y, basis){
  ncol = dim(Y)[2]
  nbas = dim(basis)[1]
  res = c()

  # Project each one separately
  for (i in 1:(ncol / nbas)){
    res = cbind(res, Y[,(i-1) * nbas + 1:nbas] %*% basis)
  }
  res
}

#' Soft impute algorithm for variables X1,X2,X3,...,Xp sparsely observed.
#' We assume they come from some continous process in \code{basis} + noise.
#'
#' @param basis orthogonal basis
#' @param K max dimension
#' @param maxIter max number of iterations
#' @param thresh threshold of for determining convergence
#' @param lambda lambda for thresholded SVD
#'
#' @return List with
#' * \code{fit} model fit. If one variable \code{fit} is a matrix. If multiple variables, \code{fit} is a list of matrices
#' * \code{u,d,v} decomposition
#' * \code{error} of the fit
#' * \code{lambda} used in the model
#'
#' @noRd
# @export
functionalMultiImpute.one = function(..., basis, K, maxIter, thresh, lambda, start = NULL){

  # arange arguments
  args <- list(...)

  Y = c()
  Yhat = c()

  for (i in 1:length(args)){
    Y = cbind(Y, args[[i]]$train)
  }

  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = 0
  err = 1e9

  if (!is.null(start) && length(args) == 1)
    Yhat = start
  if (!is.null(start) && length(args) > 1){
    Yhat = c()
    for (i in 1:length(start)){
      Yhat = cbind(Yhat, start[[i]])
    }
  }

  # Repeat SVD + impute till convergence
  dims = 1:K
  err = 0
  for (i in 1:maxIter){
    # Fill in last prediction into NULLs
    Yfill[ynas] = Yhat[ynas]

    # Run SVD on the filled matrix
    Ysvd = svd(project.on.basis(Yfill, basis))

    # Threshold SVD
    D = Ysvd$d[dims]- lambda
    D[D<0] = 0

    # Predict
    Yhat.new = project.on.basis(Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims])), t(basis))

    # Check if converged
    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-5)
    if (i %% 100 == 0)
      print(ratio)
    if (ratio < thresh)
      break
    Yhat = Yhat.new

    # Remember the error
    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
  }

  # Rearange output depending on the number of variables
  fit= list()
  ncol = dim(args[[1]]$train)[2]
  for (i in 1:length(args)){
    fit[[i]] = Yhat[,(i-1)*ncol + 1:ncol]
  }
  if (length(args) == 1){
    fit = fit[[1]]
  }

  list(fit=fit, d=D, u=Ysvd$u[,dims,drop=FALSE], id=rownames(args[[i]]), grid=as.numeric(colnames(args[[i]])), err=err, lambda = lambda)
}

# @export
functionalMultiImpute = function(..., basis = fc.basis(), K = ncol(basis), maxIter = 10e3, thresh = 10e-4, lambda = 0, final="soft"){
  args <- list(...)
  err = 1e9
  best = NULL
  bestK = K
  bestModel = NULL

  args.smpl = args
  args.smpl[[1]] = fc.sample(args.smpl[[1]])
  nargs = length(args.smpl)
  if (nargs > 1){
    for (i in 2:nargs){
      args.smpl[[i]] = apply.mask(args.smpl[[i]], args.smpl[[1]])
    }
  }
  args.smpl[["basis"]] = basis
  args.smpl[["K"]] = K
  args.smpl[["maxIter"]] = maxIter
  args.smpl[["thresh"]] = thresh

  cv.K = c()
  cv.err = c()

  for (l in lambda)
  {
    args.smpl[["lambda"]] = l
    model = do.call(functionalMultiImpute.one, args.smpl)

    err.new = 0
    for (i in 1:nargs){
      err.new = err.new + sqrt(mean((args.smpl[[i]]$test - model$fit[[i]])[args.smpl[[i]]$test.mask]**2))
    }

    cat(paste("Error with lambda=",l,"\t",err.new,"\n"))
    cv.K = c(cv.K, sum(model$d > 1e-5))
    cv.err = c(cv.err, err.new)

    if (err.new < err){
      err = err.new
      bestLambda = l
      bestK = sum(model$d > 1e-10)
      bestModel = model
    }
  }
  meta = data.frame(lambda=lambda, cv.K = cv.K, cv.err=cv.err)

  args.smpl[["lambda"]] = bestLambda
  for (i in 1:nargs){
    args.smpl[[i]]$train = args[[i]]
  }
  if (final=="hard"){
    args.smpl[["lambda"]] = 0
    args.smpl[["K"]] = bestK
    args.smpl[["start"]] = bestModel$fit
  }

  res = do.call(functionalMultiImpute.one, args.smpl)
  res$meta = meta
  res
}
