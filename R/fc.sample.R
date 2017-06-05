#' @export
fc.sample = function(X,perc=0.1){
  true = !is.na(X)
  total = sum(true)
  counts = rowSums(true)
  tosample = which(counts > 1)

  ntest = min(floor(perc * total), length(tosample))

  test = true
  test[] = FALSE
  test.rows = sample(tosample)[1:ntest]
  for (i in test.rows)
  {
    j = sample(which(true[i,]))[1]
    test[i,j] = TRUE
  }
  X.test = X
  X.train = X
  X.train[test] = NA
  X.test[!test] = NA
  list(test = X.test, train = X.train, test.rows = test.rows, test.mask = test)
}
