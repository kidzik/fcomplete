#' @export
fc.wide2long = function(X)
{
  # TODO: This function should be optimized

  time = as.numeric(colnames(X))
  ids = rownames(X)
  long = data.frame(id = c(), time = c(), value=c())

  for (i in 1:length(ids)){
    observed = X[i,!is.na(X[i,])]
    id = rep(ids[i], length(observed))
    tm = time[!is.na(X[i,])]
    toadd = data.frame(id=id,
                       time=tm,
                       value=observed)
    long = rbind(long,toadd)
  }
  rownames(long) = NULL
  long
}
