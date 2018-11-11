#' Apply a mask generated through fc.sample
#' to another matrix X of the same size
#'
#' @return list of the same format as \code{fc.sample}
#'
#' @noRd
#' @export
apply.mask = function(wide, mask){
  smp.X = list(test.mask = mask$test.mask)
  smp.X$train = wide
  smp.X$test = wide
  smp.X$test.rows = mask$test.rows
  nas = is.na(smp.X$train)
  smp.X$train[smp.X$test.mask] = NA
  smp.X$test[!smp.X$test.mask] = NA
  smp.X
}

