# @export
functionalMultiImpute = function(..., basis = fc.basis(), K = ncol(basis), maxIter = 1e4, thresh = 1e-5, lambda = 0, final="soft", mask = NULL, verbose = 1){
  args <- list(...)
  err = 1e9
  best = NULL
  bestK = K
  bestModel = NULL
  meta = NULL

  args.smpl = args

  if (!is.null(mask)){
    args.smpl[[1]] = apply.mask(args.smpl[[1]], mask)$train
  }

  args.smpl[[1]] = fc.sample(args.smpl[[1]], 0.05)
  nargs = length(args.smpl)
  if (nargs > 1){
    for (i in 2:nargs){
      args.smpl[[i]] = apply.mask(args.smpl[[i]], args.smpl[[1]])
    }
  }
  args.smpl[["basis"]] = basis
  args.smpl[["K"]] = K
  args.smpl[["maxIter"]] = maxIter
  args.smpl[["thresh"]] = thresh
  args.smpl[["verbose"]] = verbose

  cv.K = c()
  cv.err = c()
  fit.err = c()

  if (length(lambda) > 1){
    for (l in lambda)
    {
      args.smpl[["lambda"]] = l

      model <- do.call(functionalMultiImpute.one, args.smpl)

      err.new = 0
      for (i in 1:nargs){
        err.new = err.new + sqrt(mean((args.smpl[[i]]$test - model$multiFit[[i]])[args.smpl[[i]]$test.mask]**2))
      }

      print(model$err)
      if (verbose > 0)
        cat(paste("Error with lambda=",l,"\t",err.new,"\n"))
      cv.K = c(cv.K, sum(model$d > 1e-5))
      cv.err = c(cv.err, err.new)
      fit.err = c(fit.err, model$err)

      if (err.new < err){
        err = err.new
        bestLambda = l
        bestK = sum(model$d > 1e-10)
        bestModel = model
      }
    }
    meta = data.frame(lambda=lambda, cv.K = cv.K, cv.err=cv.err, fit.err = fit.err)
  }
  else {
    bestLambda = lambda
  }



  args.smpl[["lambda"]] = bestLambda
  for (i in 1:nargs){
    args.smpl[[i]]$train[args.smpl[[i]]$test.mask] = args.smpl[[i]]$test[args.smpl[[i]]$test.mask]
  }
  if (final=="hard"){
    args.smpl[["lambda"]] = 0
    args.smpl[["K"]] = bestK
    args.smpl[["start"]] = bestModel$fit
  }

  res = do.call(functionalMultiImpute.one, args.smpl)
  res$meta = meta
  res$err.cv = err
  res
}
