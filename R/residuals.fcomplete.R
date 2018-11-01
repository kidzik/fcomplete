#' @export
residuals.fcomplete = function(model){
  preds = predict(model, newdata=model$data)
  model$data[[model$y.var]] - preds
}
