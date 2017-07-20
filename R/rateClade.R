#' rateClade function
#'
#' rateClade function rescales branches according to rate (used internally)
#' @param phy_tree
#' @param branches
#' @param rate default = 1
#' @keywords 
#' @return phy_tree	
#' @export
#' @examples
#' rateClade()

rateClade <- function(phy_tree, branches, rate=1) {
	for (i in branches) phy_tree$edge.length[i] <- phy_tree$edge.length[i] * rate
	phy_tree	
}
