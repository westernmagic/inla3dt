library(INLA)
library(inlamesh3d)
library(tidyverse)
library(ggplot2)
library(assertthat)
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

# Set parameters
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

# FEM
## Spatial
spatial_fem <- mesh2fem(spatial_mesh, order = 3L)

# Model
model <- inla.spde2.matern(
	spatial_mesh,
	param = param2.matern(
		spatial_mesh,
		alpha = 2,
		prior_range = range_s_prior,
		prior_sigma = sigma_prior
	)
)

# Define a subset of locations for simulating observations
i <- sample.int(nrow(spatial_mesh$loc), ceiling(obs_ratio * nrow(spatial_mesh$loc))) %>%
	# sort for convenience (`A` in `stack` will be diagonal)
	sort()
loc <- spatial_mesh$loc[i,]

# Get precision matrix
Q <- inla.spde2.precision(
	model,
	theta = log(c(range_s_0, sigma_0) / c(range_s_prior, sigma_prior))
)

# Take a sample from the precision matrix
x <- inla.qsample(1, Q)[, 1]
invisible(assert_that(min(x) <= 0 && 0 <= max(x)))

# Projection matrix
## Spatial
A <- inla.mesh3d.make.A(spatial_mesh, loc)

# Generate noisy observations
y <- as.vector(A %*% x) + rnorm(nrow(loc), sd = sd_0)

# Stack data for INLA
stack <- inla.stack(
	data          = list(y = y),
	A             = list(A),
	effects       = list(inla.spde.make.index('field', n.spde = nrow(spatial_mesh$loc))),
	compress      = FALSE,
	remove.unused = FALSE
)

est <- inla(
	y ~ -1 + f(field, model = model),
	data              = inla.stack.data(stack),
	control.compute   = list(graph = TRUE, dic = TRUE),
	control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
	control.mode      = list(theta = log(c(sd_0^-2, range_s_0, sigma_0)) + 0.5, restart = TRUE),
	# inla.mode         = "experimental",
	control.inla      = list(int.strategy = "eb"),
	verbose           = TRUE
)
