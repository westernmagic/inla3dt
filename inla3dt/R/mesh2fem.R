#' Construct finite element method matrices.
#' 
#' Patches \code{\link[inlamesh3d]{inla.mesh3d.fem}} to support \code{order} argument
#' 
#' @param mesh
#'   an INLA mesh (\code{\link[INLA]{inla.mesh.1d()}}, \code{\link[INLA]{inla.mesh()}} or \code{\link[inlamesh3d]{inla.mesh.3d()}} object)
#' @param order
#'   The model order
#' @returns
#'   A list of sparse matrices based on basis functions \eqn{\psi_{i}}:
#'   \describe{
#'     \item{\code{c0}}{\eqn{\tilde{C}_{i, j} = \langle\psi_{i}, 1}\rangle}
#'     \item{\code{c1}}{\eqn{C_{i, j} = \langle\psi_{i}, \psi_{j}\rangle}}
#'     \item{\code{g1}}{\eqn{G^{(1)}_{i, j} = \langle\nabla\psi_{i}, \nabla\psi_{j}\rangle}}
#'     \item{\code{g}k}{\eqn{G^{(k)}_{i, j} = G^{(k - 1)}_{i, j} \tilde{C}_{i, j}^{-1} G^{(k - 1)}_{i, j}}, up to and including \code{order} = \eqn{k}}
#'   }
#'   For 1D meshes, also:
#'   \describe{
#'     \item{\code{m0}}{\eqn{M_{0} = \tilde{C}}}
#'     \item{\code{m1}}{\eqn{M_{1} =} TODO}
#'     \item{\code{m2}}{\eqn{M_{2} = G^{(1)}}}
#'   }
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' @importFrom INLA       inla.mesh.1d.fem
#' @importFrom INLA       inla.mesh.fem
#' @importFrom inlamesh3d inla.mesh3d.fem
#' @importFrom Matrix     sparseMatrix
#' @importFrom glue       glue
#'
#' @export
mesh2fem <- function(mesh, order = 2L) {
	assert_that(inherits(mesh, c("inla.mesh.1d", "inla.mesh", "inla_mesh_3d")))
	assert_that(is.scalar(order))
	assert_that(is.integer(order) || is.double(order))
	assert_that(trunc(order) == order)
	assert_that(order >= 2)
	
	if (inherits(mesh, "inla.mesh.1d")) {
		fem <- INLA::inla.mesh.1d.fem(mesh)
		# Temporal mesh matrices
		# from (Bakka et al. 2020), appendix 3.1, equations (36) and (38)
		fem$m0 <- fem$c0
		fem$m1 <- sparseMatrix(
			i = c(1, dim(fem$c0)[1]),
			j = c(1, dim(fem$c0)[1]),
			x = c(0.5, 0.5)
		)
		fem$m2 <- fem$g1
	} else if (inherits(mesh, "inla.mesh")) {
		fem <- INLA::inla.mesh.fem(mesh, order)
	} else if (inherits(mesh, "inla_mesh_3d")) {
		fem <- inlamesh3d::inla.mesh3d.fem(mesh)
		fem$c0 <- as(fem$c0, "diagonalMatrix")
		if (order >= 3) {
			for (i in 3:order) {
				fem[[glue("g{i}")]] <- fem[[glue("g{i - 1}")]] %*% (fem$g1 / fem$va)
			}
		}
	}
	
	fem
}

#' @exportS3Method inlamesh3d::inla.mesh.fem
inla.mesh.fem.inla_mesh_3d <- function(mesh, ...) {
	mesh2fem(mesh, ...)
}