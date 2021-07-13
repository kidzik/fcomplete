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

  if (length(args) > 0 && verbose > 0){
    cat(paste("Combining",length(args),"variables"))
  }
  for (i in 1:length(args)){
    Y = cbind(Y, args[[i]]$train)
  }

  yobs = !is.na(Y)
  Yres = Y
  Yres[] = 0
  Yhat = project.on.basis(Y, basis)
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
  err = sqrt(mean( ((mean(Y,na.rm = TRUE) - Y)[yobs])**2))

  for (i in 1:maxIter){
    Yhat.curves = project.on.basis(Yhat, t(basis))

    # Compute residuals of the latest prediction
    Yres[yobs] = Y[yobs] - Yhat.curves[yobs]

    # get gradient
    projected.res = project.on.basis(Yres, basis)

    # weight the gradient
    weights = sqrt(ncol(basis))/10

    # Run SVD on the current solution + gradient
    Ysvd = svd(Yhat + weights*projected.res)

    # Threshold SVD
    D = Ysvd$d[dims] - lambda
    D[D<0] = 0

    # Predict
    Yhat.new = Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims]))

    # Check if converged
    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-15)
    if (verbose && (i %% 100 == 0)){
      print(thresh)
      print(ratio)
    }
    if (ratio < thresh){
      break
    }
    Yhat = Yhat.new

  }
  numIter = i
  Yhat.curves = project.on.basis(Yhat, t(basis))

  # Remember the error
  err = sqrt( mean( (Yhat.curves - Y)[yobs]**2))

  # Rearange output depending on the number of variables
  fit= list()
  ncol = dim(args[[1]]$train)[2]
  for (i in 1:length(args)){
    fit[[i]] = Yhat.curves[,(i-1)*ncol + 1:ncol]
  }

  v = NULL
  if (ncol(t(Ysvd$v[,1:ncol(basis),drop=FALSE])) == nrow(t(basis)))
    v = t(Ysvd$v[,1:ncol(basis),drop=FALSE]) %*% t(basis)

  list(multiFit=fit, numIter = numIter, fit = fit[[1]], d=D, u=Ysvd$u[,dims,drop=FALSE], v=v[dims,,drop=FALSE], id=rownames(args[[i]]), grid=as.numeric(colnames(args[[i]])), err=err, lambda = lambda, data=args)
}
