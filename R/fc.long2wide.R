#' @export
fc.long2wide = function(groups, time, values, bins=100)
{
  # TODO: This function should be optimized

  minval = min(time)
  maxval = max(time)
  # cat(minval," ",maxval," ",(maxval - minval) / bins,"\n")
  tobin = function(x){
    ceiling ( ((x - minval) / (maxval - minval)) * bins)
  }
  binned = tobin(time)
  ids = unique(groups)
  wide = matrix(NA, length(ids), bins)

  for (i in 1:length(groups)){
    binned[i] = max(1, binned[i])
    wide[which(ids == groups[i]),binned[i]] = values[i]
  }
  rownames(wide) = ids
  colnames(wide) = seq(minval + (maxval - minval) / bins, maxval, length.out = bins)
  wide
}
