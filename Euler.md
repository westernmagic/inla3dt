```sh
env2lmod
module load \
	gcc/8.2.0 \
	git \
	r \
	cmake \
	freetype libpng zlib \
	fontconfig uuid libxml2 \
	harfbuzz glib pcre \
	fribidi \
	libtiff \
	libjpeg-turbo \
	hdf5 libszip \
	openblas

R --no-save --quiet <<-EOF
	# Set CRAN Mirror to ETHZ
	local({
		r <- getOption("repos")
		r["CRAN"] <- "https://stat.ethz.ch/CRAN"
		r["INLA"] <- "https://inla.r-inla-download.org/R/stable"
		options(repos = r)
	})

	# Install build dependencies
	install.packages(c(
		"devtools",
		"remotes",
		"Rcpp",
		"BH",
		"BiocManager"
	))

	BiocManager::install(c(
		"graph",
		"Rgraphviz"
	))
	
	install.packages(c(
		"assertthat",
		"geometry",
		"hdf5r",
		"glue",
		"dplyr",
		"xml2",
		"INLA"
	))

	remotes::install_github(
		"finnlindgren/inlamesh3d"
	)
	# Currently not working
	remote::install_github(
		"eliaskrainski/INLAspacetime"
	)

	install.packages(c(
		"tidyverse",
		"ggplot2"
	))

	# Install gmshr
	roxygen2::roxygenize(
		"gmshr",
		load_code = function(path) {
			pkgload::load_all(
				path,
				compile         = FALSE,
				helpers         = FALSE,
				attach_testthat = FALSE
			)$env
		}
	)
	Rcpp::compileAttributes("gmshr")
	devtools::document("gmshr")
	devtools::build("gmshr")
	devtools::install("gmshr")

	# Install tetgenr
	roxygen2::roxygenize(
		"tetgenr",
		load_code = function(path) {
			pkgload::load_all(
				path,
				compile         = FALSE,
				helpers         = FALSE,
				attach_testthat = FALSE
			)$env
		}
	)
	Rcpp::compileAttributes("tetgenr")
	devtools::document("tetgenr")
	devtools::build("tetgenr")
	devtools::install("tetgenr")
	
	# Install meshr
	devtools::document("meshr")
	devtools::build("meshr")
	devtools::install("meshr")

	# Install inla3dt
	devtools::document("inla3dt")
	devtools::build("inla3dt")
	devtools::install("inla3dt")
EOF
```
