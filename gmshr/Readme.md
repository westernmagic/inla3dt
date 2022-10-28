# Building
```R
roxygen2::roxygenise(
	load_code = function(path) {
		pkgload::load_all(
			path,
			compile = FALSE,
			helpers = FALSE,
			attach_testthat = FALSE
		)$env
	}
)
Rcpp::compileAttributes()
devtools::document()
devtools::build()
```

# Loading
```R
devtools::load_all(export_all = FALSE)
```