#' Replicate a matrix
#' 
#' @param x
#'   a matrix
#' @param times
#'   number of times to repeat the whole matrix
#' @param each
#'   number of times to repeat each row or column
#' @param byrow
#'   whether to repeat the columns (\code{FALSE}; the default), or the rows (\code{TRUE})
#' @returns the repeated matrix
#' 
#' @seealso \code{\link[base]{rep}}
#'
#' @importFrom assertthat assert_that
#' @export
rep.matrix <- function(x, times = 1L, each = 1L, byrow = FALSE) {
	assert_that(is.matrix(x))
	assert_that(is.integer(times))
	assert_that(!is.na(times))
	assert_that(times >= 1L)
	assert_that(is.integer(each))
	assert_that(!is.na(each))
	assert_that(each >= 1L)
	assert_that(is.logical(byrow))
	assert_that(!is.na(byrow))
	
	if (times == 1L && each == 1L) {
		return (x)
	}
	if (times != 1L && each != 1L) {
		result <- rep(
			x,
			times = times,
			each  = 1L,
			byrow = byrow
		)
		result <- rep(
			result,
			times = 1L,
			each  = each,
			byrow = byrow
		)
		return (result)
	}
	
	if (times == 1L && byrow) {
		result <- matrix(
			rep(
				as.vector(x),
				each = each
			),
			nrow = each * nrow(x),
			ncol = ncol(x)
		)
	} else if (times == 1L && !byrow) {
		result <- matrix(
			rep(
				as.vector(t(x)),
				each = each
			),
			nrow  = nrow(x),
			ncol  = each * ncol(x),
			byrow = TRUE
		)
	} else if (each == 1L && byrow) {
		result <- matrix(
			rep(
				as.vector(t(x)),
				times = times
			),
			nrow  = times * nrow(x),
			ncol  = ncol(x),
			byrow = TRUE
		)
	} else if (each == 1L && !byrow) {
		result <- matrix(
			rep(
				as.vector(x),
				times = times
			),
			nrow = nrow(x),
			ncol = times * ncol(x)
		)
	}
	
	if (nrow(result) == nrow(x)) {
		rownames(result) <- rownames(x)
	}
	if (ncol(result) == ncol(x)) {
		colnames(result) <- colnames(x)
	}
	return (result)
}