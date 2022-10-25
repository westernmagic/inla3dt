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
#' @returns A list.
#'
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#'
#' @export
interpretable2spde <- function(sigma, range_s, range_t, alpha = c(t = 1, s = 2, e = 1)) {
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
	assert_that(sigma   > 0)
	assert_that(range_s > 0)
	assert_that(range_t > 0)
	
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
	
	gamma <- list()
	gamma$s <- sqrt(8 * nu$s) * range_s^-1
	gamma$t <- range_t *(8 * (alpha$t - 0.5))^-0.5 * gamma$s^alpha$s
	gamma$e <- sigma^-1 * c_1^0.5 * gamma$t^-0.5 * gamma$s^-(alpha$alpha - 1)
	
	assert_that(gamma$t > 0)
	assert_that(gamma$s > 0)
	assert_that(gamma$e > 0)
	
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