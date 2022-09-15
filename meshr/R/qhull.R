#' Generate a mesh using qhull
#'
#' @importFrom geometry delaunayn
#' @export
qhull <- function(
	geometry,
	...
) {
	list(
		geometry = geometry,
		topology = delaunayn(geometry, ...)
	)
}