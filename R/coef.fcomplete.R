#' @export
coef.fcomplete = function(model){
  list(components = model$v, scores = model$u, weights = model$d, means = model$cmeans)
}
