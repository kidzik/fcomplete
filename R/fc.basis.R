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
fc.basis = function(d = 11, type=c("fourier","splines","linear"), norder = 4, dgrid = 100, rangeval = c(0,1)){
  if (length(type) > 1)
    type = "splines"

  if (type == "fourier")
    basis = fda::create.fourier.basis(rangeval, d)
  if (type == "splines")
    basis = fda::create.bspline.basis(rangeval, d, norder)
  if (type == "linear")
    basis = fda::create.polygonal.basis(rangeval)

  grid = seq(rangeval[1],rangeval[2],length.out=dgrid)

  S = fda::eval.basis(evalarg = grid, basisobj = basis)
  list(
    discrete = svd(S)$u, # make it orthonormal
    basis = basis
  )
}
