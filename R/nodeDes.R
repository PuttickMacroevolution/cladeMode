#' nodeDes
#'
#' nodeDes
#' @param phylo input phylogeny
#' @param node 
#' @keywords 
#' @return
#' @export
#' @examples
#' nodeDes()

nodeDes <- function(phylo, node) {
  tips <- Ntip(phylo)
  nextNodes <- node
  stopPoint <- FALSE
  storeNodes <- nextNodes
  tipTimes <- BranchingTimesFossil(phylo)
  extantTips <- tipTimes[which(tipTimes < 1e-8)]
  extantTips <-  match(names(extantTips), phylo$tip.label)
  
  while(stopPoint == FALSE) {
    startPoints <- unlist(sapply(nextNodes, function(x) which(phylo$edge[,1] == x)))
    nextNodes <- phylo$edge[startPoints, 2]
    storeNodes <- c(storeNodes, nextNodes)
    stopPoint <- all(nextNodes <= tips)
  }
  
  tipClade <- storeNodes[which(storeNodes <= Ntip(phylo))]
  tippo <-  phylo$tip.label[tipClade]
  return(tippo)
}
