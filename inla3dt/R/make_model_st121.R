#' Make DEMF(1, 2, 1) model
#' 
#' @seealso \code{\link[INLAspacetime]{st121cgeneric}}
#' 
#' @param spatial_fem
#' @param temporal_fem
#' @param range_s_prior
#' @param range_t_prior
#' @param sigma_prior
#' @returns TODO
#' 
#' @importFrom assertthat    assert_that
#' @importFrom Matrix        Diagonal
#' @importFrom INLAspacetime upperPadding
#' @importFrom INLA          inla.cgeneric.define
#' @export
make_model_st121 <- function(spatial_fem, temporal_fem, range_s_prior, range_t_prior, sigma_prior) {
	assert_that(is.list(spatial_fem))
	assert_that(is.list(temporal_fem))
	
	assert_that(is.vector(range_s_prior))
	assert_that(is.double(range_s_prior))
	assert_that(length(range_s_prior) == 2)
	
	assert_that(is.vector(range_t_prior))
	assert_that(is.double(range_t_prior))
	assert_that(length(range_t_prior) == 2)
	
	assert_that(is.vector(sigma_prior))
	assert_that(is.double(sigma_prior))
	assert_that(length(sigma_prior) == 2)
	
	n_t <- nrow(temporal_fem$g1)
	n <- n_t * nrow(spatial_fem$g1)
	m1 <- Diagonal(n_t, c(1, rep(0, n_t - 2), 1))
	lmats <- upperPadding(list(
		dc  = temporal_fem$c0 %x% spatial_fem$c0,
		dg  = temporal_fem$c0 %x% spatial_fem$g1,
		dg2 = temporal_fem$c0 %x% spatial_fem$g2,
		dg3 = temporal_fem$c0 %x% spatial_fem$g3,
		mc  = m1 %x% spatial_fem$c0,
		mg1 = m1 %x% spatial_fem$g1,
		mg2 = m1 %x% spatial_fem$g2,
		hg1 = temporal_fem$g1 %x% spatial_fem$c0,
		hg2 = temporal_fem$g1 %x% spatial_fem$g1
	))
	
	inla.cgeneric.define(
		model  = "inla_cgeneric_st121_model",
		shlib  = paste0(system.file("libs", package = "INLAspacetime"), "/cgenericModels.so"),
		n      = n,
		debug  = 0L,
		prs    = range_s_prior,
		prt    = range_t_prior,
		psigma = sigma_prior,
		ii     = lmats$graph@i + 1L,
		jj     = lmats$graph@j + 1L,
		xx     = t(lmats$xx)
	)
}