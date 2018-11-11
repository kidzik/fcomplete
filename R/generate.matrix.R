#' Helper function for bimodal positive definite matrix
#'
#' @noRd
# @export
generate.matrix = function(n, d, K)
{
  r1 = (d:1) / d #c(1,0.4,0.005,0.1 * exp(-(3:(d-1))))
  r2 = (d:1) / d #c(1.3,0.2,0.005,0.1 * exp(-(3:(d-1))))
  r1[-(1:K)] = 0
  r2[-(1:K)] = 0

  # generate covariance matrices
  V = svd(matrix(rnorm(d*d),d))$v
  D = diag(r1)
  Sigma1 = V %*% D %*% t(V)
  # V = svd(matrix(rnorm(d*d),d))$v
  # D = diag(r2)
  # Sigma2 = V %*% D %*% t(V)

  # fix the mean of one group
  # adjust another group to get overall mean equal zero
  # mm = rnorm(d)*5
  # mm = rep(0,d)

  # Generate matrices
  Ycoef = rmvnorm(n = n, sigma = Sigma1, mean = rep(0, d))
  # Ycoef2 = rmvnorm(n = n, sigma = Sigma2, mean = -mm)

  # split into two groups
  # subst = 1:n > n/2
  # Ycoef = Ycoef1
  # Ycoef[subst,] = Ycoef2[subst,]
  Ycoef
}

