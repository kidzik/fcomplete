project.on.basis = function(Y, basis){
  ncol = dim(Y)[2]
  nbas = dim(basis)[1]
  res = c()

  for (i in 1:(ncol / nbas)){
    res = cbind(res, Y[,(i-1) * nbas + 1:nbas] %*% basis)
  }
  res
}

#' @export
functionalMultiImpute.one = function(..., basis, K, maxIter, thresh, lambda){
  args <- list(...)

  Y = c()

  for (i in 1:length(args)){
    Y = cbind(Y, args[[i]]$train)
  }

  ynas = is.na(Y)
  Yfill = Y
  Yhat = Y
  Yhat[] = 0
  err = 1e9


  dims = 1:K
  err = 0
  for (i in 1:maxIter){
    Yfill[ynas] = Yhat[ynas]
    Ysvd = svd(project.on.basis(Yfill, basis))
    D = Ysvd$d[dims]- lambda
    D[D<0] = 0
    Yhat.new = project.on.basis(Ysvd$u[,dims] %*% (D * t(Ysvd$v[,dims])), t(basis))
    ratio = norm(Yhat.new - Yhat,"F") / (norm(Yhat,type = "F") + 1e-5)
    if (i %% 100 == 0)
      print(ratio)
    if (ratio < thresh)
      break
    Yhat = Yhat.new

    err = sqrt(mean( ((Yhat - Y)[!ynas])**2))
  }

  fit= list()
  ncol = dim(args[[1]]$train)[2]
  for (i in 1:length(args)){
    fit[[i]] = Yhat[,(i-1)*ncol + 1:ncol]
  }

  list(fit=fit, d=D, u=Ysvd$u[,dims,drop=FALSE], id=rownames(args[[i]]), grid=as.numeric(colnames(args[[i]])), err=err, lambda = lambda)
}

#' @export
functionalMultiImpute = function(..., basis = fc.basis(), K = ncol(basis), maxIter = 10e3, thresh = 10e-4, lambda = 0, final="soft"){
  args <- list(...)
  err = 1e9
  best = NULL
  bestK = K

  args.smpl = args
  args.smpl[[1]] = fc.sample(args.smpl[[1]])
  nargs = length(args.smpl)
  if (nargs > 1){
    for (i in 2:nargs){
      args.smpl[[i]] = fcomplete::apply.mask(args.smpl[[i]], args.smpl[[1]])
    }
  }
  args.smpl[["basis"]] = basis
  args.smpl[["K"]] = K
  args.smpl[["maxIter"]] = maxIter
  args.smpl[["thresh"]] = thresh

  for (l in lambda)
  {
    args.smpl[["lambda"]] = l
    model = do.call(functionalMultiImpute.one, args.smpl)

    err.new = 0
    for (i in 1:nargs){
      err.new = err.new + sqrt(mean((args.smpl[[i]]$test - model$fit[[i]])[args.smpl[[i]]$test.mask]**2))
    }

    cat(paste("Error with lambda=",l,"\t",err.new,"\n"))

    if (err.new < err){
      err = err.new
      bestLambda = l
      bestK = sum(model$D > 1e-10)
    }
  }
  args.smpl[["lambda"]] = bestLambda
  for (i in 1:nargs){
    args.smpl[[i]]$train = args[[i]]
  }
  if (final=="hard"){
    args.smpl[["lambda"]] = 0
    args.smpl[["K"]] = bestK
  }

  do.call(functionalMultiImpute.one, args.smpl)
}

#' @export
apply.mask = function(wide, mask){
  smp.X = list(test.mask = mask$test.mask)
  smp.X$train = wide
  smp.X$test = wide
  smp.X$test.rows = mask$test.rows
  nas = is.na(smp.X$train)
  smp.X$train[smp.X$test.mask] = NA
  smp.X$test[!smp.X$test.mask] = NA
  smp.X
}

