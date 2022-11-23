#' Convert interpretable DEMF parameters to SPDE ones.
#'
#' @seealso{Lindgren et al, 2022}
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
	assert_that(1 <= d)
	
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
	
	assert_that(sigma   > 0)
	assert_that(range_s > 0)
	assert_that(range_t > 0)
	
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
	
	gamma <- list()
	# section 3, subsection 2, subsubsection 2, equations 19-21
	gamma$s <- (8 * nu$s)^0.5 * range_s^-1
	gamma$t <- range_t * (8 * nu$t)^-0.5 * gamma$s^alpha$s
	# TODO: sigma^-1 or sigma^-0.5
	gamma$e <- c_t^0.5 * c_s^0.5 * gamma$t^-0.5 * gamma$s^-nu$s * sigma^-1
	
	# log(gamma) has to be well-defined
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
