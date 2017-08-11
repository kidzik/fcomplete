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

#' @export
# functionalMultiImpute.one.old = function(..., basis, K, maxIter, thresh, lambda){
#   args <- list(...)
#   nas = list()
#   fill = list()
#   nas.all = c()
#   hat = list()
#   hat.c = c()
#
#   for (i in 1:length(args)){
#     nas[[i]] = is.na(args[[i]])
#     nas.all = cbind(nas[[i]])
#     fill[[i]] = args[[i]]
#     hat[[i]] = fill[[i]]
#     hat[[i]][] = rnorm(prod(dim(hat[[i]])))
#     hat.c = cbind(hat.c, hat[[i]])
#   }
#
#   YX = hat.c
#
#   err = 1e9
#
#   dims = 1:K
#   for (i in 1:maxIter){
#     projected = c()
#     nc.all = list()
#
#     for (i in 1:length(args)){
#       fill[[ i ]][nas[[i]]] = hat[[ i ]][nas[[i]]]
#       pr.new = fill[[ i ]] %*% basis
#       nc.all[[i]] = ncol(pr.new)
#       projected = cbind(projected, pr.new)
#     }
#
#     YXsvd = svd(projected )
#
#     D = YXsvd$d[dims] - lambda
#     D[D<0] = 0
#     hat.new.proj = YXsvd$u[,dims] %*% (D * t(YXsvd$v[,dims]))
#
#     nc = 0
#     hat.new = c()
#     for (i in 1:length(args)){
#       hat[[i]] = hat.new.proj[,(nc+1):(nc + nc.all[[i]])] %*% t(basis)
#       hat.new = cbind(hat.new, hat[[i]])
#       nc = nc + nc.all[[i]]
#     }
#
#     ratio = norm(hat.new - hat.c, "F") / (norm(hat.c, type = "F") + 10e-10)
#     print(ratio)
#     if (ratio < thresh)
#       break
#     hat.c = hat.new
#
#     err = sqrt(mean( ((hat.c - YX)[!nas.all])**2))
#     # cat(err,sum(D > 0),"\n")
#   }
#
#   # v=(t(YXsvd$v[,dims])) %*% t(mbasis),
#   res = list(D=D, U=YXsvd$u[,dims], err=err, lambda = lambda)
#
#   ncur = 0
#   res$fit = list()
#   for (i in 1:length(args)){
#     res$fit[[i]] = hat.c[,(ncur+1):(ncur + ncol(args[[i]]))]
#     ncur = ncur + ncol(args[[i]])
#   }
#   res
# }

#' #' @export
#' functionalImpute.one = function(Y, basis, K, maxIter, thresh, lambda){
#'   ynas = is.na(Y)
#'   Yfill = Y
#'   Yhat = Y
#'   Yhat[] = 0
#'   err = 1e9
#'
#'   dims = 1:K
#'   err = 0
#'   for (i in 1:maxIter){
#'     Yfill[ynas] = Yhat[ynas]
#'     Ysvd = svd(Yfill %*% basis)
#'     D = Ysvd$d[dims]- lambda
#'     D[D<0] = 0
#'     Yhat.new = Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims])) %*% t(basis)
#'     ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 10e-10)
#'     if (ratio < thresh)
#'       break
#'     Yhat = Yhat.new
#'
#'     err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
#'     # cat(err,sum(D > 0),"\n")
#'   }
#'   list(fit=Yhat, d=diag(D), u=Ysvd$u[,dims], v=(t(Ysvd$v[,dims])) %*% t(basis), err=err, lambda = lambda)
#'
#' }
#' #' @export
#' functionalImpute = function(Y, basis = fc.basis(), K = ncol(basis), maxIter = 10e3, thresh = 10e-4, lambda = 0){
#'   err = 1e9
#'   best = NULL
#'   bestK = K
#'
#'   smpl = fc.sample(Y)
#'   for (l in lambda)
#'   {
#'     model = functionalImpute.one(smpl$train, basis, K, maxIter, thresh, l)
#'     err.new = sqrt(mean((smpl$test - model$fit)[smpl$test.mask]**2))
#'     if (err.new < err){
#'       err = err.new
#'       best = model
#'       bestLambda = l
#'       bestK = sum(model$D > 1e-10)
#'     }
#'   }
#'   functionalImpute.one(Y, basis, bestK, maxIter, thresh, bestLambda)
#' }
