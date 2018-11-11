#' @export
summary.fcomplete = function(model){
  cat("fcomplete fitted using the '",model$params$method,"' method with a formula: ",as.character(model$params$formula),"",sep="")
  cat("\n  Number of fitted subjects (N):\t",nrow(model$fit) )
  cat("\n  Number of evaluation timpoints:\t",ncol(model$fit) )
  cat("\n  Maximum dimension (K):\t\t",model$params$K )
  cat("\n  Observed to N*K ratio:\t\t", nrow(model$data)  / (model$params$K * nrow(model$fit)) )

  r.sqr = 1 - sum(residuals(model)**2) / sum((model$data[[model$y.var]])**2)
  cat("\n  Estimated r-squared:\t\t\t",r.sqr )

}
