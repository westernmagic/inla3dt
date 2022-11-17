library(INLA)
library(inlamesh3d)
library(tidyverse)
library(ggplot2)
library(assertthat)
# devtools::install_github("eliaskrainski/INLAspacetime")
library(INLAspacetime)
# devtools::install("meshr")
library(meshr)
# devtools::install("inla3dt")
library(inla3dt)

# Setup domain
## Spatial domain: a cube
spatial_domain <- rbind(
	c(0, 0, 0),
	c(0, 0, 1),
	c(0, 1, 0),
	c(0, 1, 1),
	c(1, 0, 0),
	c(1, 0, 1),
	c(1, 1, 0),
	c(1, 1, 1)
)
## Temporal domain: a line
temporal_domain <- c(0, 1)

# Set parameters
## Number of timesteps
nt <- 5
## Temporal Range
range_t_0 <- 0.5
range_t_prior <- 0.8
## Spatial Range
range_s_0 <- 0.25
range_s_prior <- 0.8
## Sigma
sigma_0 <- 0.1
sigma_prior <- 0.8
## Noise (nugget) standard deviation
sd_0 <- 0.1
## Observation / mesh point ratio
obs_ratio <- 0.8

# Parameter Assertions
invisible(assert_that(1 <= nt))
invisible(assert_that(0 <= range_t_0))
invisible(assert_that(0 < range_t_prior && range_t_prior < 1))
invisible(assert_that(0 <= range_s_0))
invisible(assert_that(0 < range_s_prior && range_s_prior < 1))
invisible(assert_that(0 <= sigma_0))
invisible(assert_that(0 < sigma_prior && sigma_prior < 1))
invisible(assert_that(0 <= sd_0))
invisible(assert_that(0 < obs_ratio && obs_ratio <= 1))

# Generate mesh
## Spatial
spatial_mesh_raw <- meshr::gmsh(
	spatial_domain,
	# Mesh.MeshSizeMax = 0.1,
	# Mesh.OptimizeThreshold = 0.3,
	# Mesh.OptimizeNetgen = 1,
	# Mesh.RecombineAll = 1,
	# Mesh.Recombine3DAll = 1
)
### Convert mesh to INLA format
spatial_mesh <- inla.mesh3d(
	spatial_mesh_raw$geometry,
	spatial_mesh_raw$topology
)
## Temporal
t <- seq(
	from       = temporal_domain[1],
	to         = temporal_domain[2],
	length.out = nt
)
temporal_mesh <- inla.mesh.1d(
	loc      = t,
	interval = temporal_domain
)

# FEM
## Spatial
spatial_fem <- mesh2fem(spatial_mesh, order = 3L)
## Temporal
temporal_fem <- mesh2fem(temporal_mesh)

# Model
model <- inla.rgeneric.define(
	model        = damf_121_rgeneric,
	spatial_fem  = spatial_fem,
	temporal_fem = temporal_fem,
	spatial_d    = 3
)

# Get precision matrix
Q <- inla.rgeneric.q(
	rmodel = model,
	cmd    = "Q",
	theta  = c(
		range_t_0,
		range_s_0,
		sigma_0
	)
)

# Q <- make_Q(
# 	spatial_fem,
# 	temporal_fem,
# 	interpretable2spde(sigma = sigma_0, range_s = range_s_0, range_t = range_t_0, d = 3)$gamma
# )

# Take a sample from the precision matrix
x <- inla.qsample(1, Q)[, 1]
invisible(assert_that(min(x) <= 0 && 0 <= max(x)))

# Define a subset of locations for simulating observations
i <- sample.int(nrow(spatial_mesh$loc), ceiling(obs_ratio * nrow(spatial_mesh$loc))) %>%
	# sort for convenience (`A` in `stack` will be diagonal)
	sort()
loc <- spatial_mesh$loc[i,]

# Projection matrix
## Spatial
spatial_A <- inla.mesh3d.make.A(spatial_mesh, loc)

## Spatio-temporal
A <- INLA::inla.spde.make.A(
		A.loc      = matrix(1, temporal_mesh$n) %x% spatial_A,
		group      = rep(temporal_mesh$loc, each = nrow(loc)),
		group.mesh = temporal_mesh
	) %>%
	zapsmall() %>%
	drop0()

# Generate noisy observations
y <- as.vector(A %*% x) + rnorm(nrow(loc) * temporal_mesh$n, sd = sd_0)

# Stack data for INLA
stack <- inla.stack(
	data          = list(y = y),
	A             = list(A),
	effects       = list(inla.spde.make.index('field', n.spde = nrow(spatial_mesh$loc), n.group = temporal_mesh$n)),
	compress      = FALSE,
	remove.unused = FALSE
)

theta_0 <- log(c(sd_0^-2, range_t_0, range_s_0, sigma_0))

est <- inla(
	y ~ -1 + f(field, model = model),
	data              = inla.stack.data(stack),
	control.compute   = list(graph = TRUE, dic = TRUE),
	control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
	control.mode      = list(theta = theta_0 + 0.5, restart = TRUE),
	# inla.mode         = "experimental",
	control.inla      = list(int.strategy = "eb"),
	verbose           = TRUE
)
