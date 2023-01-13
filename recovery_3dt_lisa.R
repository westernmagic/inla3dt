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

# Set parameters
## Number of timesteps
nt <- 5

## Temporal domain: a line
#temporal_domain <- c(0, 1)
temporal_domain <- c(0, nt)


## Temporal Range
range_t_0 <- 1.5 ##nt/2
range_t_prior <- 0.5
## Spatial Range
range_s_0 <- 0.5
range_s_prior <- 0.5
## Sigma
sigma_0 <- 1
sigma_prior <- 0.5
## Noise (nugget) standard deviation
sd_0 <- 0.1
## Observation / mesh point ratio
obs_ratio <- 1.0

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

x <- inla.qsample(1, temporal_fem$c0)[, 1]
print("summary(temporal_fem$c0)")
print(summary(x))

x <- inla.qsample(1, temporal_fem$c1)[, 1]
print("summary(temporal_fem$c1)")
print(summary(x))

x <- inla.qsample(1, temporal_fem$g1)[, 1]
print("summary(temporal_fem$g1)")
print(summary(x))


# Model
model <- make_model_121(
	spatial_mesh  = spatial_mesh,
	temporal_mesh = temporal_mesh,
	range_t_0     = range_t_0,
	range_t_prior = range_t_prior,
	range_s_0     = range_s_0,
	range_s_prior = range_s_prior,
	sigma_0       = sigma_0,
	sigma_prior   = sigma_prior,
	mode          = "cgeneric",
	##verbose       = TRUE,
	debug         = FALSE
)

### Problem: in get_Q() theta is not passed on as input
# was always using same default parameters

# Get precision matrix
#Q <- get_Q(
#	model,
#	theta = log(c(
#		range_t_0,
#		range_s_0,
#		sigma_0
#	))
#)

# use instead:
gamma_param <- interpretable2spde(sigma_0, range_s_0, range_t_0, 3)$gamma

print("gamma param")
print(gamma_param)
Q <- make_Q(spatial_fem, temporal_fem, gamma_param)

#jpeg(file="Qprior.jpeg")
#image(Q)
#dev.off()

# Take a sample from the precision matrix
x <- inla.qsample(1, Q)[, 1]
print("summary(x)")
print(summary(x))

# shouldn't be necessary, you can try without it
x = mean(x) - x

print("summary(x)")
print(summary(x))

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

jpeg(file="A.jpeg")
image(A)
dev.off()

# Generate noisy observations
noise <- rnorm(nrow(loc) * temporal_mesh$n, sd = sd_0)
y <- as.vector(A %*% x) + noise

# Stack data for INLA
# previous effects call generating incorrect matrix for A (separate space & time call is for separable models)
stack <- inla.stack(
        data          = list(y = y),
        A             = list(A),
###     effects       = list(inla.spde.make.index('field', n.spde = nrow(spatial_mesh$loc), n.group = temporal_mesh$n)),
        effects       = list(field=1:ncol(A)),
        compress      = FALSE,
        remove.unused = FALSE
)

theta_0 <- log(c(sd_0^-2, range_s_0,range_t_0, sigma_0))

est <- inla(
	y ~ -1 + f(field, model = model),
	data              = inla.stack.data(stack),
	#control.compute   = list(graph = TRUE, dic = TRUE),
	control.compute   =list(config=TRUE), keep=TRUE,
	control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
	control.mode      = list(theta = theta_0 + 0.5, restart = TRUE),
	# updated to contain noise observations not sigma_0
	control.family    = list(hyper = list(theta = list(prior = "pc.prec", param = c(sd_0^-2, 0.5)))),
	#inla.mode         = "experimental",
	control.inla      = list(int.strategy = "eb"),
	verbose           = TRUE
)

## check output, compute comparisons
gamma_0 <- interpretable2spde(sigma_0, range_s_0, range_t_0, 3)$gamma

QpriorINLA <- est$misc$configs$config[[1]]$Qprior
QpriorINLA[1:10,1:10]

cat("summary(diag(QpriorINLA))")
print(summary(diag(QpriorINLA)))

theta_INLA <- est$misc$configs$config[[1]]$theta
expThetaINLA <- exp(theta_INLA)

gamma_paramINLA <- interpretable2spde(expThetaINLA[4], expThetaINLA[2], expThetaINLA[3], 3)$gamma

names(gamma_paramINLA) <- c("t", "s", "e")
Qprior <- make_Q(spatial_fem, temporal_fem, gamma_paramINLA)

cat("summary(diag(Qprior))")
print(summary(diag(Qprior)))

Qxy_INLA <- est$misc$configs$config[[1]]$Q
Qxy_INLA[1:10,1:10]

cat("summary(diag(Qxy_INLA))")
print(summary(diag(Qxy_INLA)))

Qxy = Qprior + exp(theta_INLA[1])*t(A) %*% A

cat("summary(diag(Qxy))")
print(summary(diag(Qxy)))

Qdiff = Qxy - Qprior;
print(summary(diag(Qdiff)))

Qdiff_INLA = Qxy_INLA - QpriorINLA;
print(summary(diag(Qdiff_INLA)))


#jpeg(file="Qdiff_INLA.jpeg")
#image(Qdiff_INLA)
#dev.off()

#jpeg(file="QpriorINLA.jpeg")
#image(QpriorINLA)
#dev.off()

cat("gamma_0    :", c(exp(theta_0[1]), gamma_0), "\n")
cat("gamma INLA :", c(exp(theta_INLA[1]), gamma_paramINLA), "\n")


cat("theta_0    : ", theta_0, "\n")
cat("theta INLA : ", t(theta_INLA), "\n")

cat("exp(theta_0)    : ", theta_0, "\n")
cat("exp(theta INLA) : ", t(theta_INLA), "\n")



# Convert parameters back to interpretable scale
hyperparam <- est$summary.hyperpar[, c("0.025quant", "0.5quant", "0.975quant")] %>%
	as.data.frame() %>%
	rownames_to_column('var') %>%
	as_tibble() %>%
	mutate(
		var = case_when(
			var == 'Precision for the Gaussian observations' ~ 'sd',
			var == 'Theta1 for field'                        ~ 'theta_t',
			var == 'Theta2 for field'                        ~ 'theta_s',
			var == 'Theta3 for field'                        ~ 'theta_e',
		)
	) %>%
	rename(
		q_0.025 = `0.025quant`,
		q_0.5   = `0.5quant`,
		q_0.975 = `0.975quant`
	) %>%
	pivot_longer(-var) %>%
	pivot_wider(names_from = var, values_from = value) %>%
	rowwise() %>%
	transmute(
		quantile = name,
		sd       = sd^-0.5,
		range_t  = spde2interpretable(theta = c(t = theta_t, s = theta_s, e = theta_e), d = 3)["range_t"],
		range_s  = spde2interpretable(theta = c(t = theta_t, s = theta_s, e = theta_e), d = 3)["range_s"],
		sigma    = spde2interpretable(theta = c(t = theta_t, s = theta_s, e = theta_e), d = 3)["sigma"]
	) %>%
	ungroup() %>%
	pivot_longer(-quantile) %>%
	group_by(name) %>%
	arrange(name, value) %>%
	mutate(
		quantile = c("q_0.025", "q_0.5", "q_0.975")
	) %>%
	ungroup() %>%
	pivot_wider(names_from = quantile, values_from = value) %>%
	mutate(
		true = case_when(
			name == "sd"      ~ sd_0,
			name == "range_t" ~ range_t_0,
			name == "range_s" ~ range_s_0,
			name == "sigma"   ~ sigma_0
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
	'recovery_3dt.h5',
	'recovery_3dt.xdmf',
	geometry        = spatial_mesh$loc,
	topology        = spatial_mesh$graph$tv,
	times           = temporal_mesh$loc,
	node_attributes = list(
		observed      = matrix(rep(as.integer((1:nrow(spatial_mesh$loc)) %in% i), temporal_mesh$n), spatial_mesh$n, temporal_mesh$n),
		eta_mean      = matrix(eta$mean,      spatial_mesh$n, temporal_mesh$n),
		eta_sd        = matrix(eta$sd,        spatial_mesh$n, temporal_mesh$n),
		eta_q_0.025   = matrix(eta$q_0.025,   spatial_mesh$n, temporal_mesh$n),
		eta_q_0.5     = matrix(eta$q_0.5,     spatial_mesh$n, temporal_mesh$n),
		eta_q_0.975   = matrix(eta$q_0.975,   spatial_mesh$n, temporal_mesh$n),
		eta_d_kl      = matrix(eta$d_kl,      spatial_mesh$n, temporal_mesh$n),
		field_mean    = matrix(field$mean,    spatial_mesh$n, temporal_mesh$n),
		field_sd      = matrix(field$sd,      spatial_mesh$n, temporal_mesh$n),
		field_q_0.025 = matrix(field$q_0.025, spatial_mesh$n, temporal_mesh$n),
		field_q_0.5   = matrix(field$q_0.5,   spatial_mesh$n, temporal_mesh$n),
		field_q_0.975 = matrix(field$q_0.975, spatial_mesh$n, temporal_mesh$n),
		field_d_kl    = matrix(field$d_kl,    spatial_mesh$n, temporal_mesh$n)
	)
)
