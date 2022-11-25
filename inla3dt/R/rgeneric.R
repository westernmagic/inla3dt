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
				range_t = exp(theta[1]),
				range_s = exp(theta[2]),
				sigma   = exp(theta[3]),
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
		lambda_t <- -log(range_t_prior) / range_t_0
		lambda_s <- -log(range_s_prior) / range_s_0
		lambda_e <- -log(sigma_prior)   / sigma_0
		
		(
			  log(lambda_t) - 1         / 2 * theta[1] - lambda_t * exp(-1         / 2 * theta[1]) + log(1         / 2)
			+ log(lambda_s) - spatial_d / 2 * theta[2] - lambda_s * exp(-spatial_d / 2 * theta[2]) + log(spatial_d / 2)
			+ log(lambda_e) +                 theta[3] - lambda_e * exp(                 theta[3])
		)
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
