#' aicc weight function
#'
#' aicc weight function calculates the aicc weight of a group of models
#' @param aicc input aicc values
#' @keywords 
#' @return aicc
#' @export
#' @examples
#' fakeAICc <- c(435, 438, 434, 440)
#' aiccWeight(fakeAICc)

aiccWeight <- function(aicc) {
	delta <- aicc - min(aicc)
	allMods <- exp(-0.5 * delta)
	weights <- allMods / sum(allMods)
	return(cbind(delta, weights))
}
