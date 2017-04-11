#' aicc function
#'
#' aicc function
#' @param x likelihood
#' @param k parameters
#' @param n sample size
#' @keywords 
#' @return aicc
#' @export
#' @examples
#' aicc()

aicc <- function(x, k, n)  -2 * x + (2 * k * (n/(n - k -1)))

