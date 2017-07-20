#' subClade function
#'
#' subClade function rescales a phylogeny subclade acording to an EB model (used internally)
#' @param phy_tree
#' @param a
#' @param node
#' @param branches
#' @param times
#' @param originTime
#' @param rate default = 1
#' @keywords 
#' @return phy_tree	
#' @export
#' @examples
#' subClade()

subClade <- function(phy_tree, a, node, branches, times, originTime, rate=1) {
	phy_tree$edge.length[branches] <- phy_tree$edge.length[branches] * rate
	for (i in branches) {
   	bl <- phy_tree$edge.length[i] 
   	age <- times[i] 
   	t1 <- originTime - age	
   	t2 <- t1 + bl          	
   	phy_tree$edge.length[i] <- ((exp(a*t2)-exp(a*t1))/(a))
	}
	phy_tree	
}
