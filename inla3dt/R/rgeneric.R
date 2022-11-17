#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' @export
damf_121_rgeneric <- function(
	cmd = c(
		"graph",
		"Q",
		"mu",
		"initial",
		"log.norm.const",
		"log.prior",
		"quit"
	),
	theta = NULL
) {
	require(assertthat)
	require(inla3dt)
	
	assert_that(is.scalar(cmd))
	assert_that(is.character(cmd))
	assert_that(cmd %in% c(
		"graph",
		"Q",
		"mu",
		"initial",
		"log.norm.const",
		"log.prior",
		"quit"
	))
	
	env <- parent.env(environment())
	assert_that("spatial_fem"  %in% names(env))
	assert_that("temporal_fem" %in% names(env))
	assert_that("spatial_d"    %in% names(env))
	
	# functions
	graph <- function(theta, ...) {
		Q(theta, ...)
	}
	
	Q <- function(theta,...) {
		make_Q(
			spatial_fem,
			temporal_fem,
			interpretable2spde(
				range_t = theta[1],
				range_s = theta[2],
				sigma   = theta[3],
				d       = spatial_d
			)$gamma
		)
	}
	
	mu <- function(...) {
		numeric(0)
	}
	
	initial <- function(...) {
		c(1, 1, 1)
	}
	
	log.norm.const <- function(...) {
		numeric(0)
	}
	
	log.prior <- function(...) {
		# TODO
	}
	
	quit <- function(...) {
		invisible()
	}
	
	# theta
	if (is.null(theta) || length(theta) == 0) {
		theta <- initial()
	}
	assert_that(is.vector(theta))
	assert_that(is.double(theta))
	assert_that(length(theta) == 3)
	
	do.call(
		cmd,
		args = list(theta = theta)
	)
}
