#' Function simulates data from the paper
#'
#' @title Simulate data
#'
#' @details The function simulates data as described in the paper.
#' Generates observations \code{X1,X2} and a \code{Y = X1 + X2 + noise}
#'
#' @seealso \code{\link{fregression}}
#' @param n number of observations
#' @param d number of dimensions
#' @param K true 'low dimension'
#' @param dgrid size of the grid
#' @param clear fraction of observations to remove
#' @param noise_mag the magnitude of noise
#' @export
fsimulate = function(
  n = 100,
  d = 7,
  K = 3,
  dgrid = 51,
  clear = 0.85,
  noise.mag = 0.1
){
  # Set up a basis
  basis = fda::create.bspline.basis(c(0,1), d, 4)
  S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)

  # GENERATE DATA
  Xcoef = generate.matrix(n,d,K) * 10
  Zcoef = generate.matrix(n,d,K) * 10
  Ycoef = generate.matrix(n,d,K) * 10

  Ztrue = Zcoef %*% t(S)
  Xtrue = Xcoef %*% t(S)
  Ytrue = (Zcoef + Xcoef) %*% t(S)

  SigmaBig = genPositiveDefMat(dgrid)$Sigma
  noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise.mag
  Xnoise = Xtrue + noise
  SigmaBig = genPositiveDefMat(dgrid)$Sigma
  noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise.mag
  Znoise = Ztrue + noise
  SigmaBig = genPositiveDefMat(dgrid)$Sigma
  noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise.mag
  Ynoise = Ytrue + noise

  # Remove (1-clear)*100% of points
  nel = prod(dim(Ytrue))
  nna = ceiling(nel * clear)

  Y.wide = Ynoise
  remove.points = sample(nel)[1:nna]
  Y.wide[remove.points] = NA
  X.wide = Xnoise
  X.wide[remove.points] = NA
  Z.wide = Znoise
  Z.wide[remove.points] = NA

  keep.rows = rowSums(!is.na(Y.wide)) > 0
  X.wide = X.wide[keep.rows,]
  Z.wide = Z.wide[keep.rows,]
  Y.wide = Y.wide[keep.rows,]
  n = nrow(X.wide)

  Ytrue = Ytrue[keep.rows,]

  # TO LONG
  time = 0:(dgrid-1)/(dgrid-1)
  subj = 1:n
  data = fc.wide2long(Y.wide,time,subj,value = "Y")
  X.long = fc.wide2long(X.wide,time,subj)
  Z.long = fc.wide2long(Z.wide,time,subj)

  data$X1 = X.long$value
  data$X2 = Z.long$value
  Yobj = list(ftrue = Ytrue, fnoisy = Ynoise, fobs = Y.wide, data = data, params = list(K = K, grid = 0:(dgrid-1)/(dgrid-1)), basis = S)
}
