# @export
functionalMultiImputeCV = function(..., basis = fc.basis(), K = ncol(basis), maxIter = 10e3, thresh = 10e-4, lambda = 0, final="soft", fold = 5, cv.ratio = 0.05, verbose = verbose){
  args <- list(...)
  meta = 0

  if ( (length(lambda) == length(K)) && (length(K) == 1) ){
    bestK = K
    bestLambda = lambda
  }
  else {
    for (i in 1:fold){
      res = functionalMultiImpute(..., basis = basis, K=K, maxIter = maxIter, thresh = thresh, lambda = lambda, final = final, verbose = verbose) # TODO: cv.ratio goes here
      meta = meta + res$meta
    }
    meta = meta / fold
    bestLambda = meta[which.min(meta[,3]),1] # best lambda
    bestK = floor(meta[which.min(meta[,3]),2]) # best K
  }

  args.smpl = args
  args.smpl[["lambda"]] = bestLambda
  #  args.smpl[[1]] = fc.sample(args.smpl[[1]], cv.ratio)
  nargs = length(args)

  args.smpl[["basis"]] = basis
  args.smpl[["K"]] = K
  args.smpl[["maxIter"]] = maxIter
  args.smpl[["thresh"]] = thresh
  args.smpl[["verbose"]] = verbose

  for (i in 1:nargs){
    args.smpl[[i]] = list() # Coercing LHS to a list
    args.smpl[[i]]$train = args[[i]]
  }
  if (final=="hard"){
    args.smpl[["lambda"]] = 0
    args.smpl[["K"]] = bestK
    args.smpl[["start"]] = res$fit
  }

  res = do.call(functionalMultiImpute.one, args.smpl)
  res$meta = meta
  res
}
