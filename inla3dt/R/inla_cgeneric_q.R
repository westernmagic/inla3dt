#' @importFrom        INLA   inla
#' @importFrom        Matrix sparseMatrix
#' @importMethodsFrom Matrix t
#' @importMethodsFrom Matrix diag
inla_cgeneric_q <- function (cmodel = NULL, theta = NULL) 
{
	stopifnot(!is.null(cmodel))
	stopifnot(inherits(cmodel$f$cgeneric, "inla.cgeneric"))
	opt <- options()
	options(warn = -1)
	on.exit(options(opt))
	result <- list()
	tfile <- tempfile()
	cmodel$f$cgeneric$debug <- FALSE
	cmodel$f$cgeneric$.q <- TRUE
	cmodel$f$cgeneric$.q.file <- tfile
	try(
		inla(
			y ~ -1 + f(one, model = cmodel),
			data = data.frame(y = NA, one = 1),
			control.mode = list(theta = theta),
			silent = 2L,
			verbose = FALSE
		),
		silent = TRUE
	)
	res <- readLines(tfile)
	unlink(tfile)
	for (i in 1:length(res)) {
		res[i] <- gsub("\t", "", res[i])
		res[i] <- gsub("[ ]+", " ", res[i])
	}
	split.char <- function(i) strsplit(res[i], "[ ]+")[[1]]
	split <- function(i) as.numeric(split.char(i))
	line <- 1
	stopifnot(split.char(line)[1] == "CGENERIC_BEGIN")
	line <- line + 2
	ntheta <- split(line)[3]
	line <- line + 1
	theta <- numeric(ntheta)
	for (i in seq_len(ntheta)) {
		theta[i] <- split(line)[6]
		line <- line + 1
	}
	line <- line + 1
	result$theta <- theta
	nelm <- split(line)[6]
	line <- line + 1
	ii <- numeric(nelm)
	jj <- numeric(nelm)
	for (i in seq_len(nelm)) {
		ij <- split(line)[c(6, 9)]
		ii[i] <- ij[1]
		jj[i] <- ij[2]
		line <- line + 1
	}
	G <- sparseMatrix(i = ii + 1, j = jj + 1)
	G <- G + t(G)
	diag(G) <- 1
	result$graph <- G
	line <- line + 2
	stopifnot(split.char(line)[1] == "optimized")
	line <- line + 1
	qq <- numeric(nelm)
	for (i in seq_len(nelm)) {
		qq[i] <- split(line)[6]
		line <- line + 1
	}
	Q <- sparseMatrix(i = ii + 1, j = jj + 1, x = qq)
	dQ <- diag(Q)
	Q <- Q + t(Q)
	diag(Q) <- dQ
	result$Q <- Q
	line <- line + 1
	nelm <- split(line)[3]
	line <- line + 1
	mu <- numeric(nelm)
	for (i in seq_len(nelm)) {
		mu[i] <- split(line)[6]
		line <- line + 1
	}
	result$mu <- mu
	line <- line + 1
	result$log.prior <- split(line)[3]
	line <- line + 2
	result$log.norm.const <- split(line)[3]
	line <- line + 1
	stopifnot(split.char(line)[1] == "CGENERIC_END")
	return(result)
}
