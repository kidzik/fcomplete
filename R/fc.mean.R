#' Compute row means using observed points
#' @noRd
# @export
fc.mean = function(X){
  ms = rowMeans(X,na.rm = TRUE)
  X[] = ms
  X
}
