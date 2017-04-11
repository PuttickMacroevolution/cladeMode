#' identify nodes of a monophyletic clade
#'
#' identify nodes of a monophyletic clade
#' @param phylo
#' @param node
#' @keywords 
#' @return storeNodes
#' @export
#' @examples
#' monoClade()

monoClade <- function(phylo, node) {
  tips <- Ntip(phylo)
  nextNodes <- node
  storeNodes <- nextNodes
  stopPoint <- FALSE
  
  while(stopPoint == FALSE) {
    startPoints <- unlist(sapply(nextNodes, function(x) which(phylo$edge[,1] == x)))
    nextNodes <- phylo$edge[startPoints, 2]
    storeNodes <- c(storeNodes, nextNodes)
    stopPoint <- all(nextNodes <= tips)
  }
  return(storeNodes)
}
