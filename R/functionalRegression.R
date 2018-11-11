#' Run sparse functional regression based on latent values
#'
#' @noRd
# @export
functionalRegression = function(Y, X, basis, lambda=0, maxIter=1e5, thresh = 1e-4, K = dim(X)[2], mask = NULL, verbose=1){
  err = 1e9
  bestLambda = NULL
  bestModel = NULL
  meta = NULL

  cv.err = c()
  fit.err = c()

  args.smpl = list()
  if (!is.null(mask)){
    args.smpl[["Y"]] = apply.mask(Y,mask)
  }
  else {
    args.smpl[["Y"]] = fc.sample(Y)
  }
  args.smpl[["X"]] = X
  nargs = length(args.smpl)
  args.smpl[["basis"]] = basis
  args.smpl[["K"]] = K
  args.smpl[["maxIter"]] = maxIter
  args.smpl[["thresh"]] = thresh

  if (length(lambda) > 1){
    for (l in lambda)
    {
      args.smpl[["lambda"]] = l
      model = do.call(functionalRegression.one, args.smpl)

      # print(args.smpl[["Y"]]$test[args.smpl[["Y"]]$test.mask]**2)
      # print(args.smpl[["Y"]]$train[args.smpl[["Y"]]$test.mask]**2)
      # print(model$fit[args.smpl[["Y"]]$test.mask]**2)
      # print((args.smpl[["Y"]]$test - model$fit)[args.smpl[["Y"]]$test.mask]**2 )
      err.new = sqrt(mean((args.smpl[["Y"]]$test - model$fit)[args.smpl[["Y"]]$test.mask]**2))

      if(verbose > 0)
        cat(paste("Error with lambda=",l,"\t",err.new,"\n"))
      cv.err = c(cv.err, err.new)
      fit.err = c(fit.err, model$err)

      if (err.new < err){
        err = err.new
        bestLambda = l
        bestModel = model
      }
    }
    meta = data.frame(lambda = lambda, cv.err = cv.err, fit.err = fit.err)
  }
  else {
    bestLambda = lambda
  }

  args.smpl[["lambda"]] = bestLambda
  args.smpl[["Y"]]$train = Y
  res = do.call(functionalRegression.one, args.smpl)
  res$meta = meta
  res$err.cv = err
  res$lambda = bestLambda
  res

}
