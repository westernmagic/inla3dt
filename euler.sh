#!/bin/bash -i

set -eo pipefail
IFS=$'\n\t'

# Use new module stack
env2lmod

# Use newer GCC
module load \
    gcc/8.2.0 \
    git \
    r \
    freetype libpng zlib \
    fontconfig uuid libxml2 \
    harfbuzz glib pcre \
    fribidi \
    libtiff \
    libiconv \
    hdf5 libszip \
    cmake

# Create user library
R --no-save --quiet <<-"EOF"
    dir.create(
        path = Sys.getenv("R_LIBS_USER"),
        showWarnings = FALSE,
        recursive = TRUE
    )
EOF
# Install renv
R --no-save --quiet <<-"EOF"
    install.packages(
        "renv",
        repos = c("CRAN" = "https://stat.ethz.ch/CRAN")
    )
EOF
# Initialize renv
R --no-save --quiet <<-"EOF"
    renv::init(
        bare = TRUE,
        bioconductor = TRUE
    )
    renv::snapshot(
        repos = c(
            "CRAN" = "https://stat.ethz.ch/CRAN",
            "INLA" = "https://inla.r-inla-download.org/R/stable"
        )
    )
EOF

R --no-save --quiet <<-"EOF"
    renv::install("devtools")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("Rcpp")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("BH")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("assertthat")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("geometry")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("hdf5r")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("glue")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("dplyr")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("xml2")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("INLA")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("finnlindgren/inlamesh3d")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("westernmagic/INLAspacetime")
    renv::snapshot()
EOF
module unload libiconv
R --no-save --quiet <<-"EOF"
    renv::install("haven")
    renv::snapshot()
EOF
module load libiconv
R --no-save --quiet <<-"EOF"
    renv::install("tidyverse")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    renv::install("ggplot2")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
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
    renv::install("./gmshr")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
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
    renv::install("./tetgenr")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    devtools::document("meshr")
    renv::install("./meshr")
    renv::snapshot()
EOF
R --no-save --quiet <<-"EOF"
    devtools::document("inla3dt")
    renv::install("./inla3dt")
    renv::snapshot()
EOF

# Save loaded modules
cat > activate.sh <<-"EOF"
    env2lmod
    module load \
        gcc/8.2.0
        git \
        r \
        freetype libpng zlib \
        fontconfig uuid libxml2 \
        harfbuzz glib pcre \
        fribidi \
        libtiff \
        libiconv \
        hdf5 libszip \
        cmake
EOF

echo 'Done!'
echo 'Remember to `source activate.sh` before executing any commands!'
