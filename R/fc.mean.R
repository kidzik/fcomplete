#' Compute row means from observed points
#' @noRd
# @export
fc.mean = function(X){
  ms = rowMeans(X,na.rm = TRUE)
  X[] = ms
  X
}

