library("devtools")
library("fcomplete")
install.packages("corpcor")
library("corpcor")

dgrid = 31
d = 7
simulation = fsimulate(dgrid = dgrid,clear = 0.9, n = 100, noise.mag = 0.3, d = d, K = 2)
X = list()
X$train = simulation$fobs

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
    Ysvd = svd(fcomplete:::project.on.basis(Yfill, basis))

    # Threshold SVD
    D = Ysvd$d[dims]- lambda
    D[D<0] = 0

    # Predict
    Yhat.new = fcomplete:::project.on.basis(Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims])), t(basis))

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

system.time(functionalMultiImpute.one(X,K = 2, basis = fcomplete:::fc.basis(d = d, dgrid = dgrid), maxIter = 10000, lambda = 0, thresh = 1e-10))

