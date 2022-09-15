#' Create a mesh from a set of points.
#' 
#' @param geometry
#'   A numeric \eqn{N × 3} matrix of points (\eqn{4 ≤ N}). Required.
#' @param hull
#'   An integer \eqn{M × 3} matrix of point indices (\eqn{4 ≤ M ≤ N}) denoting the boundary surfaces (triangles). Defaults to the convex hull if the \pkg{geometry} is available.
#' @param ...
#'   Further options to pass to gmsh. See \url{https://gmsh.info/doc/texinfo/gmsh.html#Options}
#' 
#' @returns A list containing:
#'   \describe{
#'     \item{\code{$geometry}}{a numeric \eqn{N' × 3} matrix of points}
#'     \item{\code{$topology}}{an integer \eqn{M' × 4} matrix of point indices (\eqn{1 ≤ M' ≤ N'}). Each row corresponds to one tetrahedron.}
#'   }
#'
#' @references
#'   \itemize{
#'     \item C. Geuzaine and J.-F. Remacle. \emph{\href{https://gmsh.info/doc/preprints/gmsh_paper_preprint.pdf}{Gmsh: a three-dimensional finite element mesh generator with built-in pre- and post-processing facilities.}} International Journal for Numerical Methods in Engineering 79(11), pp. 1309-1331, 2009. \doi{10.1002/nme.2579}
#'     \item \href{https://gmsh.info/}{The gmsh homepage}
#'   }
#' 
#' @importFrom assertthat assert_that
#' @export
mesh <- function(
	geometry,
	hull = NULL,
	...
) {
	assert_that(!missing(geometry))
	assert_that(is.numeric(geometry))
	assert_that(is.matrix(geometry))
	assert_that(4L <= nrow(geometry))
	assert_that(ncol(geometry) == 3L)
	geometry <- apply(geometry, 2, as.numeric)
	
	if (missing(hull) || is.null(hull)) {
		if (requireNamespace("geometry", quietly = TRUE)) {
			hull <- geometry::convhulln(geometry)
		} else {
			stop("No `hull` specified, and `geometry` package not available.")
		}
	}
	assert_that(is.numeric(hull))
	assert_that(is.matrix(hull))
	assert_that(4L <= nrow(hull))
	assert_that(ncol(hull) == 3)
	assert_that(all(trunc(hull) == hull))
	assert_that(all(1L <= hull))
	assert_that(all(hull <= nrow(geometry)))
	# assert_that(!all(!apply(hull, 1, duplicated)))
	hull <- apply(hull, 2, as.integer)
	
	options <- list(...)
	if ("Mesh.FirstNodeTag" %in% names(options)) {
		warning("Overriding Mesh.FirstNodeTag to be TRUE")
	}
	options["Mesh.FirstNodeTag"] = 1L
	
	if ("Mesh.FirstElementTag" %in% names(options)) {
		warning("Overriding Mesh.FirstElementTag to be TRUE")
	}
	options["Mesh.FirstElementTag"] = 1L
	
	if ("Mesh.Renumber" %in% names(options)) {
		warning("Overriding Mesh.Renumber")
	}
	options["Mesh.Renumber"] = TRUE
	
	if (!("General.Verbosity" %in% names(options))) {
		options["General.Verbosity"] <- 0L
		if ("verbose" %in% names(options)) {
			options["General.Verbosity"] <- options["verbose"]
			options["verbose"] <- NULL
		}
	}
	
	mesh_internal(
		geometry,
		hull,
		options
	)
}