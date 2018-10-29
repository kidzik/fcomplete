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
