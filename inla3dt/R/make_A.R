#' Convert points in Cartesian coordinates to local mesh cell barycentric coordinates.
#' 
#' Based of \code{\link[inlamesh3d]{inla.mesh3d.make.A}}
#' 
#' @inheritParams cart2bary
#' @param tol
#'   Tolerance for making the result more sparse. See \code{\link[Matrix]{drop0}}
#' @returns A \eqn{(n × N)} sparse projector matrix.
#' 
#' @seealso \code{\link[INLA]{inla.spde.make.A}}
#'
#' @importFrom Matrix sparseMatrix
#' @export
make_A <- function(cartesian, geometry, topology, tol = 1e-15) {
	barycentric <- cart2bary(cartesian, geometry, topology)
	drop0(
		sparseMatrix(
			dims = c(nrow(cartesian), nrow(geometry)),
			i = rep(seq_len(nrow(barycentric$coordinates)), ncol(barycentric$coordinates)),
			j = as.vector(topology[barycentric$element, ]),
			x = as.vector(barycentric$coordinates)
		),
		tol = tol
	)
}