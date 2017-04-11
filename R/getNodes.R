#' getNodes function
#'
#' getNodes function
#' @param x likelihood
#' @param k parameters
#' @param n sample size
#' @keywords 
#' @return aicc
#' @export
#' @examples
#' getNodes()

getNodes <- function(ttr) ttr$edge[which(ttr$edge[,2] > Ntip(ttr)), 2]

