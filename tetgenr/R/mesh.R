#' Create a mesh from a set of points.
#' 
#' @param geometry
#'   A numeric \eqn{N × 3} matrix of points (\eqn{4 ≤ N}). Required.
#' @param hull
#'   An integer \eqn{M × 3} matrix of point indices (\eqn{4 ≤ M ≤ N}) denoting the boundary surfaces (triangles). Defaults to the convex hull if the \pkg{geometry} is available.
#' @param verbose
#'   The verbosity level between 0 and 3. Optional.
#' @param max_edge_radius_ratio
#'   The maximum edge-to-radius ratio. Optional.
#' @param min_dihedral_angle
#'   The minimal dihedral angle. Optional
#' @param max_tet_volume
#'   The maximum tetrahedral volume. Optional.
#' @param optimize
#'   The optimization level between 0 and 10. Optional.
#' @param flip
#'   Enable or disable edge / face flips. Optional.
#' @param smooth
#'   Enable or disable vertex smoothing. Optional.
#' @param insert_delete
#'   Enable or disable vertex insertion / deletion. Optional.
#' 
#' @returns A list containing:
#'   \describe{
#'     \item{\code{$geometry}}{a numeric \eqn{N' × 3} matrix of points}
#'     \item{\code{$topology}}{an integer \eqn{M' × 4} matrix of point indices (\eqn{1 ≤ M' ≤ N'}). Each row corresponds to one tetrahedron.}
#'   }
#'
#' @references
#'   \itemize{
#'     \item Hang Si. 2015. \emph{"TetGen, a Delaunay-Based Quality Tetrahedral Mesh Generator".} ACM Trans. on Mathematical Software. 41 (2), Article 11 (February 2015), 36 pages. \doi{10.1145/2629697}
#'     \item \href{https://wias-berlin.de/software/index.jsp?id=TetGen}{The TetGen homepage}
#'   }
#' 
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.flag
#' @importFrom assertthat is.number
#' @importFrom assertthat is.scalar
#' @export
mesh <- function(
	geometry,
	hull                  = NULL,
	verbose               = NULL,
	max_edge_radius_ratio = NULL,
	min_dihedral_angle    = NULL,
	max_tet_volume        = NULL,
	optimize              = NULL,
	flip                  = NULL,
	smooth                = NULL,
	insert_delete         = NULL
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
	assert_that(ncol(hull) == 3L)
	assert_that(all(trunc(hull) == hull))
	assert_that(all(1L <= hull))
	assert_that(all(hull <= nrow(geometry)))
	# assert_that(!all(!apply(hull, 1, duplicated)))
	hull <- apply(hull, 2, as.integer)
	
	if (!missing(verbose) && !is.null(verbose)) {
		assert_that(is.flag(verbose) || is.number(verbose))
		assert_that(trunc(verbose) == verbose)
		assert_that(0L <= verbose)
		assert_that(verbose <= 3L)
		verbose <- as.integer(verbose)
	} else {
		verbose <- 0L
	}
	
	if (!missing(max_edge_radius_ratio) && !is.null(max_edge_radius_ratio)) {
		assert_that(is.number(max_edge_radius_ratio))
		assert_that(0 <= max_edge_radius_ratio)
	} else {
		max_edge_radius_ratio <- NULL
	}
	
	if (!missing(min_dihedral_angle) && !is.null(min_dihedral_angle)) {
		assert_that(is.number(min_dihedral_angle))
		assert_that(0 <= min_dihedral_angle)
		assert_that(min_dihedral_angle <= 360)
	} else {
		min_dihedral_angle <- NULL
	}
	
	if (!missing(max_tet_volume) && !is.null(max_tet_volume)) {
		assert_that(is.number(max_tet_volume))
		assert_that(0 < max_tet_volume)
		if (is.infinite(max_tet_volume)) {
			max_tet_volume <- NULL
		}
	} else {
		max_tet_volume <- NULL
	}
	
	if (!missing(optimize) && !is.null(optimize)) {
		assert_that(is.flag(optimize) || is.number(optimize))
		assert_that(trunc(optimize) == optimize)
		assert_that(0L <= optimize)
		assert_that(optimize <= 10L)
		optimize <- as.integer(optimize)
	} else {
		optimize <- NULL
	}
	
	if (!missing(flip) && !is.null(flip)) {
		assert_that(is.flag(flip))
	} else {
		flip <- NULL
	}
	
	if (!missing(smooth) && !is.null(smooth)) {
		assert_that(is.flag(smooth))
	} else {
		smooth <- NULL
	}
	
	if (!missing(insert_delete) && !is.null(insert_delete)) {
		assert_that(is.flag(insert_delete))
	} else {
		insert_delete <- NULL
	}
	
	options <- ""
	
	if (verbose == 0L) {
		options <- paste0(options, "Q")
	} else if (verbose == 1L) {
		options <- paste0(options, "V")
	} else if (verbose == 2L) {
		options <- paste0(options, "VV")
	} else if (verbose == 3L) {
		options <- paste0(options, "VVV")
	}
	
	options <- paste0(options, "q")
	if (!is.null(max_edge_radius_ratio)) {
		options <- paste0(options, max_edge_radius_ratio)
	}
	options <- paste0(options, "/")
	if (!is.null(min_dihedral_angle)) {
		options <- paste0(options, min_dihedral_angle)
	}
	
	if (!is.null(max_tet_volume)) {
		options <- paste0(options, "a", max_tet_volume)
	}
	
	options <- paste0(options, "O")
	if (!is.null(optimize)) {
		options <- paste0(options, optimize)
	}
	options <- paste0(options, "/")
	if (!is.null(flip) || !is.null(smooth) || !is.null(insert_delete)) {
		optimization_flags <- 0L
		if (!is.null(flip) && flip) {
			optimization_flags <- optimization_flags + 1L
		}
		if (!is.null(smooth) && smooth) {
			optimization_flags <- optimization_flags + 2L
		}
		if (!is.null(insert_delete) && insert_delete) {
			optimization_flags <- optimization_flags + 4L
		}
		assert_that(is.scalar(optimization_flags))
		assert_that(is.integer(optimization_flags))
		assert_that(0L <= optimization_flags)
		assert_that(optimization_flags <= 7L)
		options <- paste0(options, optimization_flags)
	}
	
	options <- paste0(options, "C")
	
	mesh_internal(
		geometry,
		hull,
		options
	)
}