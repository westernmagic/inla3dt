#' Convert SPDE DEMF parameters to interpretable ones
#'
#' \seealso{Bakka et al. 2020, section 3}
#'
#' @param theta
#'   Optional.
#' @param gamma
#'   Optional.
#' @param alpha
#'   A vector \eqn{\alpha = (\alpha_{t}, \alpha_{s}, \alpha_{e})}
#' @importFrom assertthat assert_that
#'
#' @export
spde2interpretable <- function(theta, gamma, alpha = c(t = 1, s = 2, e = 1)) {
	# alpha
	assert_that(!is.null(alpha))
	assert_that(is.vector(alpha))
	assert_that(is.double(alpha))
	assert_that(length(alpha) == 3)
	
	if (!is.null(names(alpha))) {
		assert_that("t" %in% names(alpha))
		assert_that("s" %in% names(alpha))
		assert_that("e" %in% names(alpha))
	} else {
		names(alpha) <- c("t", "s", "e")
	}
	alpha <- as.list(alpha)
	assert_that(alpha$t >= 0)
	assert_that(alpha$s >= 0)
	assert_that(alpha$e >= 0)
	
	# theta / gamma
	assert_that(xor(!missing(theta), !missing(gamma)))
	
	## theta
	if (!missing(theta)) {
		assert_that(is.vector(theta))
		assert_that(is.double(theta))
		assert_that(length(theta) == 3)
		if (!is.null(names(theta))) {
			assert_that("t" %in% names(theta))
			assert_that("s" %in% names(theta))
			assert_that("e" %in% names(theta))
		} else {
			names(theta) <- c("t", "s", "e")
		}
		theta <- as.list(theta)
		# TODO: check
		# assert_that(theta$t > 0)
		# assert_that(theta$s > 0)
		# assert_that(theta$e > 0)
	}
	
	# gamma
	if (missing(gamma) || is.null(gamma)) {
		gamma <- sapply(theta, exp)
	}
	
	assert_that(is.vector(gamma))
	assert_that(is.double(gamma))
	assert_that(length(gamma) == 3)
	
	if (!is.null(names(gamma))) {
		assert_that("t" %in% names(gamma))
		assert_that("s" %in% names(gamma))
		assert_that("e" %in% names(gamma))
	} else {
		names(gamma) <- c("t", "s", "e")
	}
	gamma <- as.list(gamma)
	
	# transformations
	# page 17, paragraph after equation (16)
	alpha$alpha <- alpha$e + alpha$s * (alpha$t - 0.5)
	assert_that(alpha$alpha > 1)
	
	nu <- list()
	# proposition 3.1
	nu$s <- alpha$alpha - 1
	# proposition 3.2
	nu$t <- min(alpha$t - 0.5, nu$s / alpha$s)
	
	# equation (19), using factor from equation (17)
	# TODO: @Lisa check what is correct
	c_1 <- (
		(base::gamma(alpha$t - 0.5) * base::gamma(alpha$alpha - 1)) /
		(base::gamma(alpha$t)       * base::gamma(alpha$alpha) * 8 * pi^1.5)
	)
	# equation (20), resp. (17)
	sigma <- gamma$e^-1 * c_1^0.5 * gamma$t^-0.5 * gamma$s^-(alpha$alpha - 1)
	# equation (21)
	# equation (22)
	range_s <- gamma$s^-1 * sqrt(8 * nu$s)
	range_t <- gamma$t * sqrt(8 * (alpha$t - 0.5)) * gamma$s^-alpha$s
	
	c(
		sigma   = sigma,
		range_s = range_s,
		range_t = range_t
	)
}