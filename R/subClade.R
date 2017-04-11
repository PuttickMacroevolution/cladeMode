#' subClade function
#'
#' subClade function
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
   	t1 <- originTime - age	  #distance from root
   	t2 <- t1 + bl          	#length from root plus original length
   	phy_tree$edge.length[i] <- ((exp(a*t2)-exp(a*t1))/(a))
	}
	phy_tree	
}
