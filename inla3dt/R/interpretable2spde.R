#' Convert interpretable DEMF parameters to SPDE ones.
#'
#' \seealso{Bakka et al. 2020, section 3}
#'
#' @param sigma
#'   A parameter.
#' @param range_s
#'   The spatial range.
#' @param range_t
#'   The temporal range.
#' @param d
#'   The dimensionality.
#' @returns A list.
#'
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#'
#' @export
interpretable2spde <- function(sigma, range_s, range_t, d, alpha = c(t = 1, s = 2, e = 1)) {
	# d
	assert_that(is.scalar(d))
	assert_that(is.integer(d) || is.double(d))
	assert_that(trunc(d) == d)
	assert_that(d >= 1)
	
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
	
	# params
	assert_that(is.scalar(sigma))
	assert_that(is.double(sigma))
	assert_that(is.scalar(range_s))
	assert_that(is.double(range_s))
	assert_that(is.scalar(range_t))
	assert_that(is.double(range_t))
	
	# TODO: check
	# assert_that(sigma   > 0)
	# assert_that(range_s > 0)
	# assert_that(range_t > 0)
	
	# transformations
	# page 17, paragraph after equation (16)
	alpha$alpha <- alpha$e + alpha$s * (alpha$t - 0.5)
	assert_that(alpha$alpha > (d / 2))
	
	nu <- list()
	# proposition 3.1
	nu$s <- alpha$alpha - d / 2
	# proposition 3.2
	nu$t <- min(alpha$t - 0.5, nu$s / alpha$s)
	
	# equation (19), using factor from equation (17)
	c_1 <- (
		(               base::gamma(nu$t)                       * base::gamma(nu$s)       ) /
		((4 * pi)^0.5 * base::gamma(alpha$t) * (4 * pi)^(d / 2) * base::gamma(alpha$alpha))
	)
	
	gamma <- list()
	gamma$s <- sqrt(8 * nu$s) * range_s^-1
	gamma$t <- range_t * (8 * nu$t)^-0.5 * gamma$s^alpha$s
	# gamma$e <- sigma^-1 * c_1^0.5 * gamma$t^-0.5 * gamma$s^-(alpha$alpha - 1)
	# equation (19)
	gamma$e <- sigma^-0.5 * c_1^0.5 * gamma$t^-0.5 * gamma$s^-nu$s
	
	# assert_that(gamma$t > 0)
	# assert_that(gamma$s > 0)
	# assert_that(gamma$e > 0)
	
	list(
		gamma = c(
			t = gamma$t,
			s = gamma$s,
			e = gamma$e
		),
		theta = log(c(
			t = gamma$t,
			s = gamma$s,
			e = gamma$e
		))
	)
}
