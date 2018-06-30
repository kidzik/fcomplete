#' Helper function for bimodal positive definite matrix
#'
#' @noRd
# @export
generate.matrix = function(n, d)
{
  # generate covariance matrices
  V = svd(matrix(rnorm(d*d),d))$v
  D = diag(c(1,0.4,0.005,0.1 * exp(-(3:(d-1))))) * 500
  Sigma1 = V %*% D %*% t(V)
  V = svd(matrix(rnorm(d*d),d))$v
  D = diag(c(1.3,0.2,0.005,0.1 * exp(-(3:(d-1))))) * 500
  Sigma2 = V %*% D %*% t(V)

  # fix the mean of one group
  # adjust another group to get overall mean equal zero
  mm = rnorm(d)*5

  # Generate matrices
  Ycoef1 = rmvnorm(n = n, sigma = Sigma1, mean = mm)
  Ycoef2 = rmvnorm(n = n, sigma = Sigma2, mean = -2*mm)

  # split into two groups
  subst = 1:n > n/3
  Ycoef = Ycoef1
#  Ycoef[subst,] = Ycoef2[subst,]
  Ycoef
}

# @export
fc.generate.one = function (coef = 1000, basis, noise_mag)
{
  if (length(coef) == 1)
  {
    coef = generate.matrix(n, ncol(basis))
  }
  ftrue = coef %*% t(basis)
  dgrid = nrow(basis)

  SigmaBig = genPositiveDefMat(dgrid)$Sigma
  noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise_mag

  obs = ftrue + noise
  list(coef = coef, obs = obs, ftrue = ftrue)
}

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
  noise.mag = 0.2
){
  # Set up a basis
  basis = fda::create.bspline.basis(c(0,1), d, 4)
  S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis) / sqrt(dgrid)

  # GENERATE DATA
  Xcoef = generate.matrix(n,d)
  Zcoef = generate.matrix(n,d)
  Ycoef = Xcoef + Zcoef

  Ztrue = Zcoef %*% t(S)
  Xtrue = Xcoef %*% t(S)
  Ytrue = Ycoef %*% t(S)

  SigmaBig = genPositiveDefMat(dgrid)$Sigma
  noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise.mag
  Xnoise = Xtrue + noise
  noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise.mag
  Znoise = Ztrue + noise
  noise = mvrnorm(n = n, SigmaBig, mu = rep(0,dgrid)) * noise.mag
  Ynoise = Xnoise + Znoise + noise * 4

  # Remove 99% of points
  nel = prod(dim(Ytrue))
  nna = ceiling(nel * clear)

  Y.wide = Ynoise
  remove.points = sample(nel)[1:nna]
  Y.wide[remove.points] = NA
  X.wide = Xnoise
  X.wide[remove.points] = NA
  Z.wide = Znoise
  Z.wide[remove.points] = NA

  X.wide = X.wide[rowSums(!is.na(Y.wide)) > 0,]
  Z.wide = Z.wide[rowSums(!is.na(Y.wide)) > 0,]
  Y.wide = Y.wide[rowSums(!is.na(Y.wide)) > 0,]
  n = nrow(X.wide)

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
