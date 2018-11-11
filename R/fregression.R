#' Method approximates a process from sparse observations. Suppose, that for a certain subject one or multiple observations
#' are measured at some irregular timepoints. We assume that these are noise observations of some underlying process and we want
#' to approximate this process for each individual.
#'
#' For a subject \eqn{i}, we observe \eqn{Y^{i}(t),X_1^{i}(t),...,X_p^{i}(t)} at irregular subject specific \eqn{t \in t_1,...,t_p}, where \eqn{0 < t_j < T}.
#' We can bin the time interval \eqn{[0,T]} and represent each individual as a vector of fixed length with missing values.
#' Let \code{Y, X1, ..., Xp} be such matrices. Columns correspond to timepoints and rows to subjects.
#'
#' There are mulitple methods of approximating the process \code{Y}, we can:
#' * regress \code{Y} on \code{X_1,X_2,...,X_p}, we can use sparse functional regression
#' * project each subject into latent space and impute \code{Y, X_1,X_2,...,X_p} simultaniously
#' * use only information from \code{Y}, we can use functional PCA method or functional impute.
#'
#' Function \code{fregression} enables all three scenarios. Suppose \code{data} contains information in the long format, i.e. \code{data} is a
#' matrix with \eqn{p + 3} columns, where \code{data[,1]} is a \code{subjectID}, \code{data[,2]} is \code{time},  \code{data[,3]} is a value observation of \code{Y}
#' and remaining columns are covariates \code{X1, ..., Xp}. Each row corresponds to one observation for one subject.
#'
#' There are three possible \code{formula}s:
#' * \code{Y:time ~ X1 + X2 | subjectID} executes functional regression
#' * \code{1:time ~ Y + X1 + X2 | subjectID} executes dimensionality reduction
#' * \code{Y:time ~ 1 | subjectID} executes functional impute or functional PCA depending on the choice of \code{method} parameter
#'
#' @title Approximate low-rank processes from sparse longitudinal observations
#'
#' @param formula formula describing the linear relation between processes and indicating time and grouping variables. See details
#' @param data data in the long format. Use \code{\link{fc.long2wide}} and \code{\link{fc.wide2long}} for conversions
#' @param bins number of bins for matrix representation of the data
#' @param method method for functional impute: \code{fpca} for functional principal components,
#' \code{mean} for mean impute and \code{fimpute} for functional impute
#' @param lambda lambdas for SVD regularization in functional impute
#' @param d dimensionality of the basis
#' @param K upper bound of dimensionality for SVD regularization
#' @param lambda.reg lambdas for SVD regularization in regression
#' @param K.reg upper bound of dimensionality for regression
#' @param thresh thershold for convergence in functional imputee
#' @param final should the final model use \code{"hard"} or \code{"soft"} impute after choosing the optimal \code{lambda}
#' @param fold number of repetitions in cross-validation
#' @param fold how many folds in cross-validation
#' @param projection "joint" or "separate" (default). If multiple regressors are available project them jointly or separately
#' @return Returns a list
#' * \code{fit} fitted matrix \code{Y}
#' * \code{meta} results of cross-validation
#' * \code{u,d,v} svd of the underlying processes if the functional impute method has been chosen
#'
#' In case of multidimensional SVD and simultanious approximation of \code{Y,X1,X2,...,Xp}, \code{$fit} is a list of models for \code{Y,X1,X2,...,Xp}.
#' @references
#' James, Gareth M., Trevor J. Hastie, and Catherine A. Sugar.
#' \emph{Principal component models for sparse functional data.}
#' Biometrika 87.3 (2000): 587-602.
#' @examples
#' # SIMULATE DATA
#' simulation = fsimulate(seed = 1)
#' data = simulation$data
#' ftrue = simulation$ftrue
#' K = simulation$params$K
#'
# REGRESSION
#' model.mean = fregression(Y ~ time | id, data,
#'                          method = "mean")
#' model.fpca = fregression(Y ~ time | id, data,
#'                          lambda = 0, K = c(3,4,5), thresh = 1e-7, method = "fpcs")
#'
#' lambdas = c(2,3,4,5,6,8,10,12,15,20)
#' model.fimp = fregression(Y ~ time | id, data,
#'                          lambda = lambdas, thresh = 1e-5, final = "hard")
#' model.fcmp = fregression(Y + X1 + X2 ~ time | id, data, covariates,
#'                          lambda = lambdas, K = K, final = "hard")
#' model.freg = fregression(Y ~ U1 + U2 + time | id, data, model.fcmp$u,
#'                          lambda = lambdas, thresh = 1e-5,
#'                          lambda.reg = 0.1, method = "fpcs", K = K)
#' @export
fregression = function(formula, data, covariates = NULL,
                       bins = 51, method = c("fimpute", "fpcs", "mean"), lambda = c(0), maxIter = 1e5,
                       lambda.reg = 0, d = 7, K = NULL, K.reg = NULL, thresh = 1e-5, final="soft", fold = 5, cv.ratio = 0.05,
                       projection = "separate", verbose = 0, scale.covariates = TRUE)
{
  if (length(method) > 1)
    method = "fimpute"

  if (is.null(K))
    K = d
  basis = fc.basis(d = d, dgrid = bins)

  vars = parse.formula(formula)
  subj.var = vars$groups
  params = list(formula = formula, method = method, lambda = lambda, maxIter = maxIter,
                lambda.reg = lambda.reg, d = d, K = K, K.reg = K.reg, thresh = thresh, final=final, fold = fold, cv.ratio = cv.ratio,
                projection = projection)

  time.var =  vars$covariates
  time.var = time.var[time.var %in% names(data)][1]
  time = data[[time.var]]

#  if (length(vars$response) == 1)
  {
    # Functional impute or regression
    y.var = vars$response[1]

    Y = na.omit(data[,c(subj.var, time.var, y.var)])
    Y.wide = fc.long2wide(Y[,1], as.numeric(Y[,2]), as.numeric(Y[,3]), bins = bins)
    params[["Y.wide"]] = Y.wide

    # Estimate population mean
    LE = lowess(time, data[[y.var]])
    cmeans = approx(x = LE$x, y = LE$y, xout = seq(min(time),max(time),length.out = bins))$y

    Y.wide = t(t(Y.wide) - cmeans)
#    yscale = sd(Y.wide[!is.na(Y.wide)])
    Y.wide[!is.na(Y.wide)] = Y.wide[!is.na(Y.wide)] #/ yscale
  }

  mint = min(time)
  maxt = max(time)
  time.grid = ((0:(bins - 1)) / (bins-1)) * (maxt - mint) + mint


  # Case 1: Y ~ time : id -- do functional impute
  if (length(vars$response) == 1 && length(vars$covariates) == 1)
  {
    if (method == "fpcs"){
      res = fc.fpca(Y, d = d, K = K, grid.l = 0:(bins-1)/(bins-1))
      # res$fit = t(t(res$fit) - cmeans) / yscale # silly but consistent
    }
    else if (method == "mean"){
      res = list(fit = fc.mean(Y.wide), v = 0)
    }
    else {
      res = functionalMultiImputeCV(Y.wide, basis = basis, lambda = lambda, K = K, thresh = thresh, final = final, fold = fold, cv.ratio = cv.ratio, maxIter = maxIter, verbose = verbose)
    }
    if (method != "fpcs"){
      res$fit = t(t(res$fit) + cmeans)
    }
    res$Y = t(t(Y.wide) + cmeans)
    res$cmeans = cmeans
    res$basis = basis

    res$time.grid = time.grid
    res$subj.var = subj.var
    res$y.var = y.var
    res$time.var = time.var
    res$data = data
    res$params = params

    if (method != "mean"){
      rownames(res$v) = paste("eigenfunction",1:nrow(res$v))
    }
    row.names(res$fit) = row.names(Y.wide)
    class(res) = "fcomplete"
    return(res)
  }

  # Assume there are covariates
  X.long = list()
  X.wide = list()

  nvars = length(vars$response)
  for (i in 1:nvars)
  {
    X.long[[i]] = na.omit(data[,c(subj.var, time.var, vars$response[i])])
    X.wide[[i]] = fc.long2wide(X.long[[i]][,1], as.numeric(X.long[[i]][,2]), as.numeric(X.long[[i]][,3]), bins = bins)

    if (scale.covariates){
      X.long[[i]][,3] = scale(X.long[[i]][,3], scale = FALSE)
      X.wide[[i]][!is.na(X.wide[[i]])] = scale(X.wide[[i]][!is.na(X.wide[[i]])], scale = FALSE)
    }
  }
  params[["X.wide"]] = X.wide

  # Case 2: X + Y ~ time -- do unsupervised learning
  if (length(vars$response) > 1 && length(vars$covariates) == 1)
  {
    args = X.wide
    args$basis = basis
    args$lambda = lambda
    args$thresh = thresh
    args$final = final
    args$K = K
    args$verbose = verbose
    res = do.call(functionalMultiImputeCV, args)
    row.names(res$fit) = row.names(Y.wide)

    res$time.grid = time.grid
    res$subj.var = subj.var
    res$y.var = y.var
    res$time.var = time.var
    res$data = data
    res$params = params
    colnames(res$u) = paste("U",1:ncol(res$u),sep="")
    rownames(res$u) = row.names(Y.wide)

    class(res) = "fcomplete"
    return(res)
  }

#   # Case 3: Y ~ Y + X -- do principal component regression with Y
#   models = list()
#   combinedU = c()
#
#
#   for (i in 1:nvars)
#   {
#     # # skip response, it will be used separately
#     # if (length(vars$response) == 2 && vars$response[1] != vars$covariates[i])
#     # {
#     X.long[[i]][,3] = scale(X.long[[i]][,3], scale = FALSE)
#     X.wide[[i]][!is.na(X.wide[[i]])] = scale(X.wide[[i]][!is.na(X.wide[[i]])], scale = FALSE)
#     if (method == "fpcs"){
#       models[[i]] = fc.fpca(X.long[[i]], d = d, K = K, grid.l = 0:(bins-1)/(bins-1))
#       combinedU = cbind(combinedU, models[[i]]$fpcs)
#     }
#     else {
# #     models[[i]] = functionalMultiImpute(X.wide[[i]], basis = basis, lambda = lambda, thresh = thresh, K = K, final = final, mask = maskedY)
# #     models[[i]] = functionalMultiImputeCV(X.wide[[i]], basis = basis, lambda = lambda, K = K, thresh = 0, final = final, fold = fold, cv.ratio = cv.ratio, maxIter = maxIter)
# #     print(names(X.long[[i]]))
#       nm = names(X.long[[i]])
#      models[[i]] = fregression(paste0(nm[3],":",nm[2]," ~ 1|",nm[1]), X.long[[i]], bins = bins, lambda = lambda, K = K, thresh = thresh, final = final, fold = fold, cv.ratio = cv.ratio, maxIter = maxIter, method="fimpute", verbose = verbose)
#      combinedU = cbind(combinedU, models[[i]]$u)
#     }
#     # }
#     # else {
#     #   print(paste("Skipping",vars$response[1]))
#     # }
#   }
#   if (method == "fimpute" && projection == "joint"){
#     args.smpl = X.wide
#     args.smpl[["basis"]] = basis
#     args.smpl[["lambda"]] = lambda
#     args.smpl[["thresh"]] = thresh
#     args.smpl[["K"]] = K
#     args.smpl[["final"]] = final
#     res = do.call(functionalMultiImpute, args.smpl)
#     combinedU = res$u
#   }
#
#   if (is.null(K.reg))
#     K.reg = ncol(Y.wide)

#  maskedY = fc.sample(Y.wide, 0.05)
  combinedU = covariates[,colnames(covariates) %in% vars$covariates]
  combinedU = cbind(1,scale(combinedU))

  # Case 3: Y ~ X -- do principal component regression without Y
  res = functionalRegression(Y.wide, combinedU, basis, lambda = lambda, K = K, thresh = 1e-10, maxIter = maxIter)
  res$Y = t(t(Y.wide) + cmeans)
  res$X = X.wide
  res$U = combinedU
#  res$X.models = models
  res$fit = t(t(res$fit) + cmeans)

  res$time.grid = time.grid
  res$subj.var = subj.var
  res$y.var = y.var
  res$time.var = time.var
  res$data = data
  res$params = params

  row.names(res$fit) = row.names(Y.wide)
  class(res) = "fcomplete"
  return(res)
}

# case4 = function(){
#   # Case 4 experimental: Y ~ Y + X -- do principal component regression with Y
#   Y.tmp = Y.wide
#
#   lastFitR = 0
#   lastFitI = 0
#
#   lastfit = 1
#   for (i in 1:200){
#     resR = functionalRegression(Y.tmp, combinedU, basis, K = 1, thresh = 1e-10, mask = maskedY, verbose = 0)
#     resR$fit = resR$fit*0.3
#     Y.tmp = Y.wide - resR$fit
#     resI = functionalMultiImpute(Y.tmp, basis = basis, thresh = thresh, verbose = 0, lambda = lambda.reg[1], K = 1) #, final = final, fold = fold, cv.ratio = cv.ratio, maxIter = maxIter)
#     #    resI = fc.fpca(Y.tmp, d=d, K = 2, grid.l = 0:(bins-1)/(bins-1))
#     Y.tmp = Y.wide - resI$fit
#
#     dR = norm(resR$fit - lastFitR, type = "F") / norm(resR$fit, type="F")
#     dI = norm(resI$fit - lastFitI, type = "F") / norm(resI$fit, type="F")
#
#     residuum = Y.wide - resI$fit
#     residuum[is.na(residuum)] = 0
#     print(norm(residuum, type='F'))
#
#     residuum = Y.wide - resR$fit
#     residuum[is.na(residuum)] = 0
#     print(norm(residuum, type='F'))
#
#     residuum = Y.wide - resR$fit - resI$fit
#     residuum[is.na(residuum)] = 0
#     res.norm = norm(residuum, type='F')
#     print(res.norm)
#
#     stopcond = abs(res.norm - lastfit) / lastfit
#     print(stopcond)
#
#     if (stopcond < 0.0001)
#       break
#     lastfit = res.norm
#
#     #    if (dR + dI < 0.2)
#     #      break
#
#     lastFitR = resR$fit
#     lastFitI = resI$fit
#   }
#   res = list()
#   sm = resR$fit + resI$fit
#   res$cmeans = cmeans
#   res$fit = t(t(sm) + cmeans)
#   res$fitI = resI$fit
#   res$fitR = resR$fit
#   res$U = combinedU
#   res$reps = i
#   res
# }
