#' Rescale phylogeny with early burst values
#'
#' Internal function to rescale phylogeny with early burst values
#' @param phy_tree input phylogeny
#' @param a 
#' @param node
#' @param branches
#' @param times
#' @param originTime 
#' @keywords 
#' @return rescaled vcv matrix
#' @export
#' @examples
#' cladeMode()


ebSubClade <- function (phy_tree, a, node, branches, times, originTime) {
  for (i in branches) {
    bl <- phy_tree$edge.length[i]
    age <- times[i]
    t1 <- originTime - age	  #distance from root
    t2 <- t1 + bl          	#length from root plus original length
    phy_tree$edge.length[i] <- (exp(a*t2)-exp(a*t1))/(a)
  }
  split.vcv.node(phy_tree, node)[[2]]
}