#' @export
fitted.fcomplete = function(model){
  predict(model, newdata=model$data)
}
