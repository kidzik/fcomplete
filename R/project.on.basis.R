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
