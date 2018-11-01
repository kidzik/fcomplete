#' @export
summary.fcomplete = function(model){
  cat("fcomplete fitted using the",model$params$method,"method" )
  cat("\nNumber of fitted subjects:\t",nrow(model$fit) )
  cat("\nNumber of evaluation timpoints:\t",ncol(model$fit) )
}
