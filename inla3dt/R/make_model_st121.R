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
make_model_st121 <- function(
	spatial_fem,
	temporal_fem,
	range_s_prior,
	range_t_prior,
	sigma_prior,
	d,
	alpha = c(t = 1, s = 2, e = 1),
	debug = FALSE
) {
	assert_that(is.list(spatial_fem))
	assert_that(is.list(temporal_fem))
	
	# params
	assert_that(is.vector(range_s_prior))
	assert_that(is.double(range_s_prior))
	assert_that(length(range_s_prior) == 2)
	
	assert_that(is.vector(range_t_prior))
	assert_that(is.double(range_t_prior))
	assert_that(length(range_t_prior) == 2)
	
	assert_that(is.vector(sigma_prior))
	assert_that(is.double(sigma_prior))
	assert_that(length(sigma_prior) == 2)
	
	# d
	assert_that(is.scalar(d))
	assert_that(is.integer(d) || is.double(d))
	assert_that(trunc(d) == d)
	
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
	
	assert_that(
		(alpha$t == 1 && alpha$s == 0 && alpha$e == 2) ||
		(alpha$t == 1 && alpha$s == 2 && alpha$e == 1) ||
		(alpha$t == 2 && alpha$s == 0 && alpha$e == 2) ||
		(alpha$t == 2 && alpha$s == 2 && alpha$e == 0)
	)
	
	# debug
	assert_that(is.scalar(debug))
	assert_that(is.logical(debug))
	
	
	# transformations
	# page 17, paragraph after equation (16)
	alpha$alpha <- alpha$e + alpha$s * (alpha$t - 0.5)
	assert_that(alpha$alpha > 1)
	
	nu <- list()
	# proposition 3.1
	nu$s <- alpha$alpha - 1
	# proposition 3.2
	nu$t <- min(alpha$t - 0.5, nu$s / alpha$s)
	
	# alphas <- as.integer(strsplit(model, '')[[1]])
	# nu.t <- alphas[1]-1/2
	# alpha <- alphas[3] + alphas[2]*nu.t
	# nu.s <- alpha-1
	
	cc <- log(c(
		(8 * nu$s)^0.5,
		(8 * nu$t)^-0.5,
		(
			(               base::gamma(alpha$t - 0.5)                    * base::gamma(alpha$alpha - 1)) /
			((4 * pi)^0.5 * base::gamma(alpha$t)       * (4 * pi)^(d / 2) * base::gamma(alpha$alpha)    )
		)	
		# 0.5*(lgamma(nu.t) - lgamma(alphas[1]) -1.5*log(4*pi))
	))
	
	mm <- stModel.matrices(smesh, tmesh, model)
	n <- smesh$n * tmesh$n
	nm <- ncol(mm$TT)
	stopifnot(nm == length(mm$bb))
	jmm <- pmatch(paste0('M', 1:nm), names(mm))
	stopifnot(length(jmm[complete.cases(jmm)]) == nm)
	
	llib <- system.file("libs", package = "INLAspacetime")
	lmats <- upperPadding(mm[jmm], relative = FALSE)
	stopifnot(n == nrow(lmats$graph))
	
	return(do.call(
		"inla.cgeneric.define",
		list(
			model    = "inla_cgeneric_sstspde",
			shlib    = paste0(llib, "/INLAspacetime.so"),
			n        = n,
			debug    = 0L,
			ii       = lmats$graph@i,
			jj       = lmats$graph@j,
			aaa      = alphas,
			manifold = as.integer(manifold),
			nm       = as.integer(nm),
			cc       = as.double(cc),
			bb       = mm$bb,
			prs      = range_s_prior,
			prt      = range_t_prior,
			psigma   = sigma_prior,
			tt       = t(mm$TT),
			xx       = t(lmats$xx)))
		)
	
	# n_t <- nrow(temporal_fem$g1)
	# n <- n_t * nrow(spatial_fem$g1)
	# m1 <- Diagonal(n_t, c(1, rep(0, n_t - 2), 1))
	# lmats <- upperPadding(list(
	# 	dc  = temporal_fem$c0 %x% spatial_fem$c0,
	# 	dg  = temporal_fem$c0 %x% spatial_fem$g1,
	# 	dg2 = temporal_fem$c0 %x% spatial_fem$g2,
	# 	dg3 = temporal_fem$c0 %x% spatial_fem$g3,
	# 	mc  = m1 %x% spatial_fem$c0,
	# 	mg1 = m1 %x% spatial_fem$g1,
	# 	mg2 = m1 %x% spatial_fem$g2,
	# 	hg1 = temporal_fem$g1 %x% spatial_fem$c0,
	# 	hg2 = temporal_fem$g1 %x% spatial_fem$g1
	# ))
	# 
	# inla.cgeneric.define(
	# 	model  = "inla_cgeneric_sstspde_model",
	# 	shlib  = paste0(system.file("libs", package = "INLAspacetime"), "/INLAspacetime.so"),
	# 	n      = n,
	# 	debug  = ifelse(debug, 1L, 0L),
	# 	prs    = range_s_prior,
	# 	prt    = range_t_prior,
	# 	psigma = sigma_prior,
	# 	ii     = lmats$graph@i + 1L,
	# 	jj     = lmats$graph@j + 1L,
	# 	xx     = t(lmats$xx)
	# )
}
