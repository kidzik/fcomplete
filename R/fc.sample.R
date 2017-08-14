# @export
fc.sample = function(X, perc = 0.1, one.per.row=TRUE){
  true = !is.na(X)
  total = sum(true)
  counts = rowSums(true)
  tosample = which(counts > 1)

  if (one.per.row){
    ntest = min(floor(perc * total), length(tosample))

    test = true
    test[] = FALSE
    test.rows = sample(tosample)[1:ntest]
    for (i in test.rows)
    {
      j = sample(which(true[i,]))[1]
      test[i,j] = TRUE
    }
  }
  else {
    nel = prod(dim(true))
    nna = ceiling(nel * perc)
    test = true
    test[] = FALSE
    remove.points = sample(nel)[1:nna]
    test[remove.points] = TRUE
    test.rows = which(rowSums(test) > 0)
  }

  X.test = X
  X.train = X
  X.train[test] = NA
  X.test[!test] = NA
  list(test = X.test, train = X.train, test.rows = test.rows, test.mask = test)
}
