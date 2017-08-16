#' Creates an orthonormal basis
#'
#' @param d dimension
#' @param type fourier or splines
#' @param norder order of splines
#' @param dgrid size of the grid
#'
#' @return matrix with a basis orthonormalized by svd
#' @noRd
# @export
fc.basis = function(d = 11, type=c("fourier","splines"), norder = 4, dgrid = 100){
  if (length(type) > 1)
    type = "splines"

  if (type == "fourier")
    basis = fda::create.fourier.basis(c(0,1), d)
  if (type == "splines")
    basis = fda::create.bspline.basis(c(0,1), d, norder)

  S = fda::eval.basis(evalarg = 0:(dgrid-1)/(dgrid-1), basisobj = basis)
  svd(S)$u # make it orthonormal
}
