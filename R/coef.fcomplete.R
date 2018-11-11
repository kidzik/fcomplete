#' @export
coef.fcomplete = function(model){
  regression = NULL
  if("res.reg" %in% names(model))
    regression = model$res.reg$coef %*% t(model$basis)
  list(components = model$v, scores = model$u, weights = model$d, means = model$cmeans, regression = regression )
}
