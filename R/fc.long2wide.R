# @export
fc.long2wide = function(groups, time, values, bins=100, ids = NULL, time.lim = c(min(time),max(time)))
{
  # TODO: This function should be optimized

  ord = order(groups, time)
  groups = groups[ord]
  time = time[ord]
  values = values[ord]

  if (is.null(ids)){
    ids = unique(groups)
  }

  # scale the time and bin the timepoints
  minval = min(time)
  maxval = max(time)
  tobin = function(x){
    ceiling ( ((x - minval) / (maxval - minval)) * bins)
  }
  binned = tobin(time)

  # create the output matrix of the right dimensions
  wide = matrix(NA, length(ids), bins)

  # go through all observations and fill in points
  for (i in 1:length(groups)){
    binned[i] = max(1, binned[i])
    wide[which(ids == groups[i]),binned[i]] = values[i]
  }

  # annotate rows and columns
  rownames(wide) = ids
  colnames(wide) = seq(minval + (maxval - minval) / bins, maxval, length.out = bins)
  wide
}

#' @export
sample.long = function(X, id.var, time.var, value.var, ratio = 0.1, min.per.sbj = 2){
#  X = X[sample(1:nrow(X)),] # shuffle
#  X = X[order(X[[id.var]]),]  # sort by id (keep the rest shuffled)
  res = list(X = X)
  ntest = nrow(X)*ratio
  tbl.X = table(X[[id.var]])
  nm = names(tbl.X[tbl.X > min.per.sbj])
  test.mask = X[[id.var]] %in% nm
  test.obs = sample(which(test.mask))
  test.obs = test.obs[1:min(ntest,length(test.obs))]
  test.mask[] = FALSE
  test.mask[test.obs] = TRUE

  res$test.mask = test.mask
  res$test.obs = test.obs
  res$train = res$X[!res$test.mask,]
  res$test = res$X[res$test.mask,]
  res$test.matrix = fc.long2wide(groups = res$test[[id.var]], time = res$test[[time.var]],
                                 values = res$test[[value.var]], bins = 51,
                                 ids = unique(X[[id.var]]), time.lim = c(min(res$X[[time.var]]),max(res$X[[time.var]])))
  res$train.matrix = fc.long2wide(groups = res$train[[id.var]], time = res$train[[time.var]],
                                 values = res$train[[value.var]], bins = 51,
                                 ids = unique(X[[id.var]]), time.lim = c(min(res$X[[time.var]]),max(res$X[[time.var]])))
  res$train.matrix = res$train.matrix[rowSums(!is.na(res$train.matrix)) > 0.5,]
  res$test.matrix = res$test.matrix[row.names(res$train.matrix),]
  res
}
