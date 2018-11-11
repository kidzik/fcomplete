#' For a given matrix X, split it to a train and test set.
#' If \code{one.per.row} is \code{TRUE}, sample
#' at most 1 point from each row
#'
#' @noRd
# @export
fc.sample = function(X, perc = 0.1, one.per.row=TRUE)
{
  true = !is.na(X)

  if (one.per.row){
    # choose rows to sample from
    total = sum(true)
    counts = rowSums(true)
    tosample = which(counts > 1)

    # generate the test mask by iterating through rows with
    # at least two observations
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
    # take all the observations and sample the test mask
    nel = prod(dim(true))
    nna = ceiling(nel * perc)
    test = true
    test[] = FALSE
    remove.points = sample(nel)[1:nna]
    test[remove.points] = TRUE
    test.rows = which(rowSums(test) > 0)
  }

  # create the test set from the sampled mask
  # and mask the train set
  X.test = X
  X.train = X
  X.train[test] = NA
  X.test[!test] = NA

  list(test = X.test, train = X.train, test.rows = test.rows, test.mask = test)
}
