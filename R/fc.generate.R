# Helper function for bimodal positive definite matrix
generate.matrix = function(n, d)
{
  V = svd(matrix(rnorm(d*d),d))$v
  D = diag(c(1,0.9,0.5,exp(-(3:(d-1))))) * 500
  Sigma1 = V %*% D %*% t(V)
  V = svd(matrix(rnorm(d*d),d))$v
  D = diag(c(1.3,0.4,0.4,exp(-(3:(d-1))))) * 500
  Sigma2 = V %*% D %*% t(V)
  Ycoef1 = rmvnorm(n = n, sigma = Sigma1, mean = rnorm(d))
  Ycoef2 = rmvnorm(n = n, sigma = Sigma2, mean = rnorm(d))
  subst = runif(n) > 0
  Ycoef = Ycoef1
  Ycoef[subst,] = Ycoef2[subst,]
  Ycoef
}

#' @export
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
