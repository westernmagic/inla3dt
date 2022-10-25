#' Make precision matrix
#' 
#' @param spatial_fem
#'   The spatial FEM
#' @param temporal_fem
#'   The temporal FEM
#' @param gamma
#'   The SPDE parameters
#' @returns
#'   the precision matrix
#' 
#' @seealso Bakka 2020, appendix 3.2, equation (56)
#'
#' @importFrom assertthat assert_that
#' @export
make_Q <- function(spatial_fem, temporal_fem, gamma) {
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
	
	gamma$e^2 * (
		temporal_fem$c0 %x% make_Q_s(spatial_fem, gamma$s, 3)
		+ 2 * gamma$t   * temporal_fem$m1 %x% make_Q_s(spatial_fem, gamma$s, 2)
		+     gamma$t^2 * temporal_fem$g1 %x% make_Q_s(spatial_fem, gamma$s, 1)
	)
}