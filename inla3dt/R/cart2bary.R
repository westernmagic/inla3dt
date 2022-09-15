#' Convert points in Cartesian coordinates to local mesh cell barycentric coordinates.
#' 
#' Based of \code{\link[inlamesh3d]{inla.mesh3d.bary}}
#' 
#' @param cartesian
#'   A \eqn{(n × d)}-\link[base]{matrix} of points in \eqn{d}-dimensional Cartesian space to be converted.
#' @param geometry
#'   A \eqn{(N × d)}-\link[base]{matrix} of mesh vertices in \eqn{d}-dimensional Cartesian coordinates.
#' @param topology
#'   A \eqn{(M × (d + 1))}-\link[base]{matrix} of mesh cell indices.
#' @returns A list containing
#'   \describe{
#'     \item{\code{$coordinates}}{the barycentric coordinates}
#'     \item{\code{$elements}}{to which element the barycentric coordinates refer to}
#'   }
#' 
#' @seealso \code{\link[geometry]{cart2bary}}
#'
#' @importFrom assertthat assert_that
cart2bary <- function(cartesian, geometry, topology) {
	assert_that(is.matrix(cartesian))
	assert_that(is.matrix(geometry))
	assert_that(is.matrix(topology))
	assert_that(is.double(cartesian) || is.integer(cartesian))
	assert_that(is.double(geometry) || is.integer(geometry))
	assert_that(is.integer(topology))
	
	assert_that(ncol(cartesian) == ncol(geometry))
	d <- ncol(cartesian)
	
	assert_that(ncol(topology) == d + 1)
	assert_that(1 <= min(topology))
	assert_that(max(topology) <= nrow(geometry))
	n <- nrow(cartesian)
	
	# pre-allocate memory
	cartesian1 <- cbind(cartesian, 1)
	barycentric <- matrix(-Inf, n, d + 1)
	elements <- integer(n)
	
	# for each element
	for (e in seq_len(nrow(topology))) {
		# get element vertices
		vertices <- geometry[topology[e, ], , drop = FALSE]
		vertices1 <- cbind(vertices, 1)
		# solve least squares
		w <- solve(t(vertices1), t(cartesian1))
		# check which points lie closer to the element than previously
		i <- which(apply(w, 2, min) > apply(barycentric, 1, min))
		if (length(i) > 0) {
			barycentric[i,] <- t(w[, i, drop = FALSE])
			elements[i] <- e
		}
	}
	
	list(coordinates = barycentric, elements = elements)
}