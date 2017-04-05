find.lambda = function(Psi, S){
  Psi.inv = sqrt(solve(Psi))

  S.svd = svd(Psi.inv %*% S %*% Psi.inv)

  V = (S.svd$v)
  D = S.svd$d - 1
  D[D < 0] = 0
  sqrt(Psi) %*% V %*% diag(sqrt(D))
}

mm.fit = function(Y,B){
  X = Y %*% B

  sigma = 1
  S = cov(X)

  ncol = dim(B)[2]
  to.optim = function(sigma) {
    Psi = sigma**2 * diag(ncol)
    Lambda = find.lambda(Psi, S)
    res = sum(tr((Lambda) %*% t(Lambda) + Psi - S)**2)
    print(paste(sigma,res))
    res
  }

  res = optim(1, to.optim, method = "Brent",lower = 0.00001, upper = max(abs(S)))
  list(Lambda = find.lambda(res$par**2 * diag(ncol), S), sigma = res$par)
}

mm.scores = function(M, Y, B){
  ncol = dim(M$Lambda)[1]
  Psi = diag(rep(M$sigma ** 2, ncol))
  t(solve(diag(ncol) +  t(M$Lambda) %*% solve(Psi) %*% M$Lambda) %*% t(M$Lambda) %*% solve(Psi) %*% t(Y %*% B))
}
