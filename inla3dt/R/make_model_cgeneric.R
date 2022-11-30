#' Make DEMF(1, 2, 1) model
#' 
#' @seealso \code{\link[INLAspacetime]{st121cgeneric}}
#' 
#' @param spatial_mesh
#'   the spatial mesh
#' @param temporal_mesh
#'   the temporal mesh
#' @param range_t_0
#' @param range_t_prior
#' @param range_s_0
#' @param range_s_prior
#' @param sigma_0
#' @param sigma_prior
#' @param mode
#' @param ...
#' @returns TODO
#' 
#' @importFrom assertthat    assert_that
#' @export
make_model_121 <- function(
	spatial_mesh,
	temporal_mesh,
	range_t_0,
	range_t_prior,
	range_s_0,
	range_s_prior,
	sigma_0,
	sigma_prior,
	mode = "cgeneric",
	...
) {
	assert_that(is.scalar(mode))
	assert_that(is.character(mode))
	assert_that(mode %in% c("cgeneric", "rgeneric"))
	
	if (mode == "cgeneric") {
		stModel.define(
			smesh          = spatial_mesh,
			tmesh          = temporal_mesh,
			model          = "121",
			control.priors = list(
				prt    = c(range_t_0, range_t_prior),
				prs    = c(range_s_0, range_s_prior),
				psigma = c(sigma_0,   sigma_prior)
			),
			useINLAprecomp = FALSE,
			...
		)
	} else if (mode == "rgeneric") {
		inla.rgeneric.define(
			model         = demf_121_rgeneric,
			spatial_fem   = spatial_fem,
			temporal_fem  = temporal_fem,
			spatial_d     = 3,
			range_t_0     = range_t_0,
			range_t_prior = range_t_prior,
			range_s_0     = range_s_0,
			range_s_prior = range_s_prior,
			sigma_0       = sigma_0,
			sigma_prior   = sigma_prior,
			...
		)
	}
}
