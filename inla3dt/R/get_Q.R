#' Get precision matrix from model
#' 
#' @param model
#'   the model
#' @param theta
#'   theta for rgeneric
#' @return the precision matrix
#' 
#' @importFrom assertthat assert_that
#' @export
get_Q <- function(
	model,
	theta = NULL
) {
	assert_that(inherits(model, c("inla.cgeneric", "inla.rgeneric")))
	
	if (inherits(model, "inla.cgeneric")) {
		inla.cgeneric.q(model)$Q
	} else if (inherits(model, "inla.rgeneric")) {
		inla.rgeneric.q(model, "Q", theta)
	}
}