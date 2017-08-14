# @export
fc.mean = function(X){
  ms = rowMeans(X,na.rm = TRUE)
  X.true = X
#  X = t(X)
  X[] = ms
#  X = t(X)
#  X[!is.na(X.true)] = X.true[!is.na(X.true)]
  X
}
