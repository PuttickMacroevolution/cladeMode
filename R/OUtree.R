#' rateClade function
#'
#' rateClade function (internal function) rescale tree according to an OU process using alpha
#' @param phy_tree
#' @param alpha
#' @param branches
#' @param times
#' @param originTime
#' @keywords 
#' @return phy_tree	
#' @export
#' @examples
#' OUtree()

OUtree <- function (phy_tree, alpha, branches, times, originTime) {
    Tmax <- originTime
    for (i in branches) {
       	bl <- phy_tree$edge.length[i]
   	age <- times[i]
	t1 = originTime - age
	t2 = t1 + bl
	phy_tree$edge.length[i] <- (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t2)) * (1 - exp(-2 * alpha * t2)) - (1/(2 * alpha)) * exp(-2 * alpha * (Tmax - t1)) * (1 - exp(-2 * alpha * t1)) }
	return(phy_tree)
}	
