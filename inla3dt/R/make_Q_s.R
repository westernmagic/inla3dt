#' Make spatial precision matrix
#'
#' \deqn{Q_{k, i} = \sum_{j = 0}^{i} \binom{j}{i} \kappa^{2 (i - j)} G_{j}}
#'
#' @param fem
#'   The FEM
#' @param kappa
#'   The \eqn{\kappa}
#' @param order
#'   The order
#' @returns
#'   The precision matrix
#'
#' @seealso Lindgren 2011, section 2.3, equation (10)
#'
#' @importFrom assertthat assert_that
#' @importFrom assertthat is.scalar
#' @importFrom glue       glue
#' @importFrom purrr      reduce
#' @export
make_Q_s <- function(fem, kappa, order) {
	assert_that(!is.null(kappa))
	assert_that(is.scalar(kappa))
	assert_that(is.double(kappa))
	
	assert_that(!is.null(order))
	assert_that(is.scalar(order))
	assert_that(is.double(order) || is.integer(order))
	assert_that(trunc(order) == order)
	assert_that(order >= 1)
	
	assert_that(!is.null(fem))
	assert_that(is.list(fem))
	assert_that("c0" %in% names(fem))
	for (i in 1:order) {
		assert_that(glue("g{i}") %in% names(fem))
	}
	fem$g0 <- fem$c0
	
	reduce(sapply(0:order, function(i) {
		choose(order, i) * kappa^(2 * (i - order)) * fem[[glue("g{i}")]]
	}), `+`)
}