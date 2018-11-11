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
