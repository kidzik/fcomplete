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
functionalMultiImpute.one = function(..., basis, K, maxIter, thresh, lambda, start = NULL, verbose = 1){

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
  K = min(K,ncol(basis))
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
    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-15)
    if (verbose && (i %% 100 == 0))
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

  v = NULL
  if (ncol(t(Ysvd$v[,1:ncol(basis),drop=FALSE])) == nrow(t(basis)))
    v = t(Ysvd$v[,1:ncol(basis),drop=FALSE]) %*% t(basis)

  list(multiFit=fit, fit = fit[[1]], d=D, u=Ysvd$u[,dims,drop=FALSE], v=v, id=rownames(args[[i]]), grid=as.numeric(colnames(args[[i]])), err=err, lambda = lambda, data=args)
}

# @export
functionalMultiImpute = function(..., basis = fc.basis(), K = ncol(basis), maxIter = 1e4, thresh = 1e-5, lambda = 0, final="soft", mask = NULL, verbose = 1){
  args <- list(...)
  err = 1e9
  best = NULL
  bestK = K
  bestModel = NULL
  meta = NULL

  args.smpl = args

  if (!is.null(mask)){
    args.smpl[[1]] = apply.mask(args.smpl[[1]], mask)$train
  }

  args.smpl[[1]] = fc.sample(args.smpl[[1]], 0.05)
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
  args.smpl[["verbose"]] = verbose

  cv.K = c()
  cv.err = c()
  fit.err = c()

  if (length(lambda) > 1){
    for (l in lambda)
    {
      args.smpl[["lambda"]] = l
      model = do.call(functionalMultiImpute.one, args.smpl)

      err.new = 0
      for (i in 1:nargs){
        err.new = err.new + sqrt(mean((args.smpl[[i]]$test - model$multiFit[[i]])[args.smpl[[i]]$test.mask]**2))
      }

      print(model$err)
      if (verbose > 0)
        cat(paste("Error with lambda=",l,"\t",err.new,"\n"))
      cv.K = c(cv.K, sum(model$d > 1e-5))
      cv.err = c(cv.err, err.new)
      fit.err = c(fit.err, model$err)

      if (err.new < err){
        err = err.new
        bestLambda = l
        bestK = sum(model$d > 1e-10)
        bestModel = model
      }
    }
    meta = data.frame(lambda=lambda, cv.K = cv.K, cv.err=cv.err, fit.err = fit.err)
  }
  else {
    bestLambda = lambda
  }



  args.smpl[["lambda"]] = bestLambda
  for (i in 1:nargs){
    args.smpl[[i]]$train[args.smpl[[i]]$test.mask] = args.smpl[[i]]$test[args.smpl[[i]]$test.mask]
  }
  if (final=="hard"){
    args.smpl[["lambda"]] = 0
    args.smpl[["K"]] = bestK
    args.smpl[["start"]] = bestModel$fit
  }

  res = do.call(functionalMultiImpute.one, args.smpl)
  res$meta = meta
  res$err.cv = err
  res
}

# @export
functionalMultiImputeCV = function(..., basis = fc.basis(), K = ncol(basis), maxIter = 10e3, thresh = 10e-4, lambda = 0, final="soft", fold = 5, cv.ratio = 0.05){
  args <- list(...)
  meta = 0

  for (i in 1:fold){
    res = functionalMultiImpute(..., basis = basis, K=K, maxIter = maxIter, thresh = thresh, lambda = lambda, final = final)
    meta = meta + res$meta
  }
  meta = meta / fold

  args.smpl = args
  args.smpl[["lambda"]] = meta[which.min(meta[,3]),1]
  args.smpl[[1]] = fc.sample(args.smpl[[1]], cv.ratio)
  nargs = length(args)

  args.smpl[["basis"]] = basis
  args.smpl[["K"]] = K
  args.smpl[["maxIter"]] = maxIter
  args.smpl[["thresh"]] = thresh

  print(nargs)
  for (i in 1:nargs){
    args.smpl[[i]]$train = args[[i]]
  }
  if (final=="hard"){
    args.smpl[["lambda"]] = 0
    args.smpl[["K"]] = bestK
    args.smpl[["start"]] = res$fit
  }

  res = do.call(functionalMultiImpute.one, args.smpl)
  res$meta = meta
  res
}
