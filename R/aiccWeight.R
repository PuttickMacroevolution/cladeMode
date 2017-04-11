#' aicc weight function
#'
#' aicc weight function
#' @param aicc input aicc values
#' @keywords 
#' @return aicc
#' @export
#' @examples
#' aiccWeight()

aiccWeight <- function(aicc) {
	delta <- aicc - min(aicc)
	allMods <- exp(-0.5 * delta)
	weights <- allMods / sum(allMods)
	return(cbind(delta, weights))
}
