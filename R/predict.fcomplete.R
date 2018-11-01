predict.fcomplete.one = function(model, newdata){
  timepoint = newdata[[model$time.var]]
  subj.id = newdata[[model$subj.var]]
  # row = model$fit[rownames(model$fit) %in% subj.id,]
  #
  # n = length(model$time.grid)
  # if (max(timepoint) > model$time.grid[n] || min(timepoint) < model$time.grid[1]){
  #   error("prediction outside the time scale")
  # }
  # m = sum(model$time.grid <= timepoint)
  # row[m]
  predict.slice(model, subj.id, timepoint)
}

predict.slice = function(model, ids, time){
  subj.id = ids
  row = model$fit[rownames(model$fit) %in% subj.id,]

  n = length(model$time.grid)
  epsilon = 1e-5
  if (max(time)-epsilon > model$time.grid[n] || min(time)+epsilon < model$time.grid[1]){
    stop("prediction outside the time scale")
  }
  time.idx = c()
  for (i in 1:length(time)){
    time.idx = c(time.idx, sum(model$time.grid <= time[i] + epsilon))
  }
  res = model$fit[rownames(model$fit) %in% ids, time.idx]
  res = as.matrix(res)
  colnames(res) = time
  res
}

#' @export
predict.fcomplete = function(model, newdata = NULL, ids = NULL, time = NULL){
  # If no slicing return the entire fit
  if (is.null(ids) && is.null(newdata) && is.null(time))
    return(model$fit)

  # If newdata is given then slice by ids which are there
  if (!is.null(newdata))
    ids = unique(newdata[,model$subj.var])

  # Predict one by one
  if (!is.null(newdata)){
    pred = rep(NA, nrow(newdata))
    for (i in 1:nrow(newdata)){
      pred[i] = predict.fcomplete.one(model, newdata[i,])
    }
    #  newdata[[model$y.val]] = pred
    return(pred)
  }

  # If sliced by id but not time then return selected ids
  if (!is.null(ids) && is.null(time)){
    res = model$fit[rownames(model$fit) %in% ids,]
    colnames(res) = model$time.grid
    return(res)
  }

  # If sliced by id and time, interpolate
  if (!is.null(ids) && !is.null(time)){
    return(predict.slice(model,ids,time))
  }
}
