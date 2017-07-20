#' monoTips
#'
#' monoTips calculates the proportion of the tree diversity represented by a clade
#' @param phylo input phylogeny
#' @param node 
#' @keywords 
#' @return
#' @export
#' @examples
#' monoTips()

monoTips <- function(phylo, node) {
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
  tipCladeExtant <- length(which(complete.cases(match(extantTips, tipClade))))
  return(c(length(tipClade),  tipCladeExtant , tips - length(tipClade)))
}
