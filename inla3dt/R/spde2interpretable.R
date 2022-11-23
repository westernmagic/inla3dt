#' Convert SPDE DEMF parameters to interpretable ones
#'
#' \seealso{Bakka et al. 2020, section 3}
#'
#' @param theta
#'   Optional.
#' @param gamma
#'   Optional.
#' @param d
#'   The dimensionality
#' @param alpha
#'   A vector \eqn{\alpha = (\alpha_{t}, \alpha_{s}, \alpha_{e})}
#' @importFrom assertthat assert_that
#'
#' @export
spde2interpretable <- function(theta, gamma, d, alpha = c(t = 1, s = 2, e = 1)) {
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
	assert_that(gamma$t > 0)
	assert_that(gamma$s > 0)
	assert_that(gamma$e > 0)
	
	# section 3, subsection 2, proposition 1, page 9
	alpha$alpha <- alpha$e + alpha$s * (alpha$t - 0.5)
	assert_that(alpha$alpha > (d / 2))
	
	nu <- list()
	# section 3, subsection 2, proposition 1, page 9
	nu$s <- alpha$alpha - d / 2
	# section 3, subsection 2, proposition 2, page 9
	nu$t <- min(alpha$t - 0.5, nu$s / alpha$s)
	
	# section 3, subsection 2, equation 16
	c_s <- base::gamma(nu$s) / (base::gamma(alpha$alpha) * (4 * pi)^(d / 2))
	c_t <- base::gamma(nu$t) / (base::gamma(alpha$t)     * (4 * pi)^0.5    )
	
	# section 3, subsection 2, subsubsection 2, equation 19
	# equation (20), resp. (17)
	sigma <- c_t^0.5 * c_s^0.5 * gamma$t^-0.5 * gamma$e^-1 * gamma$s^-nu$s
	range_s <- gamma$s^-1 * (8 * nu$s)^0.5
	range_t <- gamma$t * gamma$s^-alpha$s * (8 * nu$t)^0.5
	
	c(
		sigma   = sigma,
		range_s = range_s,
		range_t = range_t
	)
}