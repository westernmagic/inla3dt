library(INLA)
library(inlamesh3d)
library(tidyverse)
library(ggplot2)
# devtools::install_github("eliaskrainski/INLAspacetime")
library(INLAspacetime)
# devtools::install_github("statguy/SpaceTimeModels")
devtools::document("meshr")
devtools::load_all("meshr", export_all = FALSE)
devtools::document("inla3dt")
devtools::load_all("inla3dt", export_all = FALSE)

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
nt <- 3
## Temporal Range
range_t_0 <- 0.5
range_t_prior <- 0.8
## Spatial Range
range_s_0 <- 0.5
range_s_prior <- 0.8
## Sigma
sigma_0 <- 4
sigma_prior <- 0.8
## Noise (nugget) standard deviation
sd_0 <- 0.1

# Generate mesh
## Spatial
spatial_mesh_raw <- meshr::gmsh(spatial_domain)
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

# Define a subset of locations for simulating observations
i <- sample.int(nrow(spatial_mesh$loc), 100) %>%
	# sort for convenience (`A` in `stack` will be diagonal)
	sort()
loc <- spatial_mesh$loc[i,]

loc <- spatial_mesh$loc

# Get precision matrix
Q <- make_Q(
	spatial_fem,
	temporal_fem,
	interpretable2spde(sigma = sigma_0, range_s = range_s_0, range_t = range_t_0)$gamma
)

# Take a sample from the precision matrix
x <- inla.qsample(1, Q)[, 1]

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
	effects       = list(inla.spde.make.index('field', n.spde = nrow(loc), n.group = temporal_mesh$n)),
	compress      = FALSE,
	remove.unused = FALSE
)

model <- make_model_st121(
	spatial_fem,
	temporal_fem,
	c(range_s_0, range_s_prior),
	c(range_t_0, range_t_prior),
	c(sigma_0, sigma_prior)
)

est <- inla(
	y ~ -1 + f(field, model = model),
	data              = inla.stack.data(stack),
	control.compute   = list(graph = TRUE, dic = TRUE),
	control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
	control.mode      = list(theta = log(c(sd_0^-2, range_s_0, range_t_0, sigma_0)) + 0.5, restart = TRUE),
	inla.mode         = "experimental",
	control.inla      = list(int.strategy = "eb"),
	verbose           = TRUE
)

# Convert parameters back to interpretable scale
est$summary.hyperpar[, c("0.025quant", "0.5quant", "0.975quant")] %>%
	as.data.frame() %>%
	rownames_to_column('var') %>%
	as_tibble() %>%
	mutate(
		var = case_when(
			var == 'Precision for the Gaussian observations' ~ 'sd',
			var == 'Theta1 for field'                        ~ 'range',
			var == 'Theta2 for field'                        ~ 'sigma'
		)
	) %>%
	rename(
		q_0.025 = `0.025quant`,
		q_0.5   = `0.5quant`,
		q_0.975 = `0.975quant`
	) %>%
	pivot_longer(-var) %>%
	mutate(
		value = case_when(
			var == 'range' ~ exp(value) * range_prior,
			var == 'sigma' ~ exp(value) * sigma_prior,
			var == 'sd'    ~ 1 / sqrt(value)
		)
	) %>%
	mutate(
		name = case_when(
			var == 'sd' & name == 'q_0.025' ~ 'q_0.975',
			var == 'sd' & name == 'q_0.975' ~ 'q_0.025',
			TRUE                            ~ name
		)
	) %>%
	arrange(name) %>%
	pivot_wider(var) %>%
	mutate(
		true = case_when(
			var == 'range' ~ range_0,
			var == 'sigma' ~ sigma_0,
			var == 'sd'    ~ sd_0
		)
	) %>%
	mutate(
		recovered = q_0.025 <= true & true <= q_0.975
	)

eta <- est$summary.linear.predictor %>%
	rownames_to_column("var") %>%
	as_tibble() %>%
	mutate(
		kind = case_when(
			str_detect(var, "^APredictor\\.") ~ "A",
			str_detect(var, "^Predictor\\.")  ~ "meshpoint"
		) %>% as.factor(),
		index = var %>% str_remove("^\\w+\\.") %>% as.integer()
	) %>%
	select(
		kind,
		index,
		mean,
		sd,
		q_0.025 = `0.025quant`,
		q_0.5   = `0.5quant`,
		q_0.975 = `0.975quant`,
		d_kl    = kld
	) %>%
	filter(kind == "meshpoint") %>%
	arrange(index)

field <- est$summary.random$field %>%
	as_tibble() %>%
	select(
		index = ID,
		mean,
		sd,
		q_0.025 = `0.025quant`,
		q_0.5   = `0.5quant`,
		q_0.975 = `0.975quant`,
		d_kl    = kld
	) %>%
	mutate(index = as.integer(index))

write_xdmf(
	'recovery.h5',
	'recovery.xdmf',
	geometry      = spatial_mesh$loc,
	topology      = spatial_mesh$graph$tv,
	observed      = as.integer((1:nrow(spatial_mesh$loc)) %in% i),
	eta_mean      = eta$mean,
	eta_sd        = eta$sd,
	eta_q_0.025   = eta$q_0.025,
	eta_q_0.5     = eta$q_0.5,
	eta_q_0.975   = eta$q_0.975,
	eta_d_kl      = eta$d_kl,
	field_mean    = field$mean,
	field_sd      = field$sd,
	field_q_0.025 = field$q_0.025,
	field_q_0.5   = field$q_0.5,
	field_q_0.975 = field$q_0.975,
	field_d_kl    = field$d_kl
)

