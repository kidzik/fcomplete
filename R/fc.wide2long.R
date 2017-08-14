# @export
fc.wide2long = function(X, time = NULL, ids = NULL, value = "value")
{
  # TODO: This function should be optimized

  if (is.null(time))
    time = as.numeric(colnames(X))
  if (is.null(ids))
    ids = rownames(X)
  long = data.frame(id = c(), time = c(), value=c())

  for (i in 1:length(ids)){
    observed = X[i,!is.na(X[i,])]
    if (length(observed) == 0){
#      stop(paste("Empty row",ids[i]))
    }
    id = rep(ids[i], length(observed))
    tm = time[!is.na(X[i,])]
    toadd = data.frame(id=id,
                       time=tm,
                       value=observed)
    long = rbind(long,toadd)
  }
  rownames(long) = NULL
  colnames(long)[3] = value
  long
}
