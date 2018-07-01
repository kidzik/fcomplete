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
#' @param K upper bound of dimensionality for SVD regularization
#' @param lambda.reg lambdas for SVD regularization in regression
#' @param K.reg upper bound of dimensionality for regression
#' @param thresh thershold for convergence in functional imputee
#' @param final should the final model use \code{"hard"} or \code{"soft"} impute after choosing the optimal \code{lambda}
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
#' model.mean = fregression(Y:time ~ 1 | id, data,
#'                          method = "mean")
#' model.fpca = fregression(Y:time ~ 1 | id, data,
#'                          lambda = 0, K = c(3,4,5), thresh = 1e-7, method = "fpcs")
#'
#' lambdas = c(2,3,4,5,6,8,10,12,15,20)
#' model.fimp = fregression(Y:time ~ 1 | id, data,
#'                          lambda = lambdas, thresh = 1e-5, final = "hard")
#' model.fcmp = fregression(0:time ~ Y + X1 + X2 | id, data,
#'                          lambda = lambdas, K = K, final = "hard")
#' model.freg = fregression(Y:time ~ X1 + X2 | id, data,
#'                          lambda = lambdas, thresh = 1e-5,
#'                          lambda.reg = 0.1, method = "fpcs", K = K)
#' @export
fregression = function(formula, data,
                       bins = 51, method = c("fimpute", "fpcs", "mean"), lambda = c(0), maxIter = 1e5,
                       lambda.reg = 0, K = NULL, K.reg = NULL, thresh = 1e-5, final="soft", d = 7, fold = 5, cv.ratio = 0.05)
{
  if (length(method) > 1)
    method = "fimpute"

  if (is.null(K))
    K = d
  basis = fc.basis(d = d, dgrid = bins)

  vars = parse.formula(formula)
  subj.var = vars$groups

  if (length(vars$response) == 1)
  {
    # If no Y and we do unsupervised learning
    time.var = vars$response[1]
  }
  else
  {
    # If there is Y to regress
    y.var = vars$response[1]
    time.var = vars$response[2]
    Y = na.omit(data[,c(subj.var, time.var, y.var)])
    Y.wide = fc.long2wide(Y[,1], as.numeric(Y[,2]), as.numeric(Y[,3]), bins = bins)

    time = data[[time.var]]
    LE = lowess(time, data[[y.var]])
    cmeans = approx(x = LE$x, y = LE$y, xout = seq(min(time),max(time),length.out = bins))$y
    Y.wide = t(t(Y.wide) - cmeans)
    yscale = sd(Y.wide[!is.na(Y.wide)])
    Y.wide[!is.na(Y.wide)] = Y.wide[!is.na(Y.wide)] / yscale
  }

  # Case 1: Y ~ 1 -- do functional impute
  if (length(vars$covariates) == 0)
  {
    if (method == "fpcs"){
      res = fc.fpca(Y, d = d, K = K, grid.l = 0:(bins-1)/(bins-1))
      res$fit = t(t(res$fit) - cmeans) / yscale # silly but consistent
    }
    else if (method == "mean"){
      res = list(fit = fc.mean(Y.wide))
    }
    else {
      res = functionalMultiImputeCV(Y.wide, basis = basis, lambda = lambda, K = K, thresh = thresh, final = final, fold = fold, cv.ratio = cv.ratio, maxIter = maxIter)
    }
    res$fit = t(t(res$fit * yscale) + cmeans)
    res$Y = t(t(Y.wide * yscale) + cmeans)
    return(res)
  }

  # Assume there are covariates
  X.long = list()
  X.wide = list()

  nvars = length(vars$covariates)
  for (i in 1:nvars)
  {
    X.long[[i]] = na.omit(data[,c(subj.var, time.var, vars$covariates[i])])
    X.wide[[i]] = fc.long2wide(X.long[[i]][,1], as.numeric(X.long[[i]][,2]), as.numeric(X.long[[i]][,3]), bins = bins)
  }

  # Case 2: 1 ~ X -- do unsupervised learning
  if (length(vars$response) == 1)
  {
    args = X.wide
    args$basis = basis
    args$lambda = lambda
    args$thresh = thresh
    args$final = final
    args$K = K
    return(do.call(functionalMultiImputeCV, args))
  }

  # Case 3: Y ~ Y + X -- do principal component regression with Y
  # Case 4: Y ~ X -- do principal component regression without X
  models = list()
  combinedU = c()

  maskedY = fc.sample(Y.wide, 0.05)

  for (i in 1:nvars)
  {
    # skip response, it will be used separately
    if (length(vars$response) == 2 && vars$response[1] != vars$covariates[i])
    {
      X.long[[i]][,3] = scale(X.long[[i]][,3])
      X.wide[[i]][!is.na(X.wide[[i]])] = scale(X.wide[[i]][!is.na(X.wide[[i]])])
      if (method == "fpcs"){
        models[[i]] = fc.fpca(X.long[[i]], d = d, K = K, grid.l = 0:(bins-1)/(bins-1))
        combinedU = cbind(combinedU, models[[i]]$fpcs)
      }
      else {
        models[[i]] = functionalMultiImpute(X.wide[[i]], basis = basis, lambda = lambda, thresh = thresh, K = K, final = final, mask = maskedY)
        combinedU = cbind(combinedU, models[[i]]$u)
      }
    }
    else {
      print(paste("Skipping",vars$response[1]))
    }
  }

  if (is.null(K.reg))
    K.reg = ncol(Y.wide)

  # Case 3: Y ~ X -- do principal component regression without Y
  if (length(vars$response) == 2 && !(vars$response[1] %in% vars$covariates)){
    print("Case 3")
    res = functionalRegression(Y.wide, combinedU, basis, lambda = lambda.reg, K = K.reg, thresh = 1e-10, mask = maskedY)
    res$Y = t(t(Y.wide * yscale) + cmeans)
    res$X = X.wide
    res$U = combinedU
    res$X.models = models
    res$fit = t(t(res$fit * yscale) + cmeans)
    return(res)
  }

  # Case 4 experimental: Y ~ Y + X -- do principal component regression with Y
  print("Case 4")
  Y.tmp = Y.wide

  lastFitR = 0
  lastFitI = 0
  for (i in 1:20){
    resI = functionalMultiImpute(Y.tmp, basis = basis, K = 1, thresh = thresh, verbose = 0) #, final = final, fold = fold, cv.ratio = cv.ratio, maxIter = maxIter)
    Y.tmp = Y.wide - resI$fit
    resR = functionalRegression(Y.tmp, combinedU, basis, K = 1, thresh = 1e-10, mask = maskedY, verbose = 0)
    Y.tmp = Y.wide - resR$fit

    dR = norm(resR$fit - lastFitR, type = "F") / norm(resR$fit, type="F")
    dI = norm(resI$fit - lastFitI, type = "F") / norm(resI$fit, type="F")
    print(c(dR, dI))

    if (dR + dI < 0.01)
      break

    lastFitR = resR$fit
    lastFitI = resI$fit
  }
  res = list()
  sm = (resI$fit + resR$fit)/2
  res$fitI = t(t(resI$fit * yscale) + cmeans)
  res$fitR = t(t(resR$fit * yscale) + cmeans)
  res$yscale = yscale
  res$cmeans = cmeans
  res$fit = t(t(sm * yscale) + cmeans)
  res
}
