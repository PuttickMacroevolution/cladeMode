#' nodeDes
#'
#' identifies the tip descendents of an ancestral node
#' @param phylo input tree in APE formate
#' @param node 
#' @keywords 
#' @return a vector of tip labels
#' @export
#' @examples
#' @examples
#' library(geiger)
#' set.seed(30)
#' # simulate tree with 20 species
#' simTree <- rcoal(20)
#' nodeDes(simTree, 22)

nodeDes <- function (phylo, node) 
{
    tips <- Ntip(phylo)
    nextNodes <- node
    stopPoint <- FALSE
    storeNodes <- nextNodes
    while (stopPoint == FALSE) {
        startPoints <- unlist(sapply(nextNodes, function(x) which(phylo$edge[, 
            1] == x)))
        nextNodes <- phylo$edge[startPoints, 2]
        storeNodes <- c(storeNodes, nextNodes)
        stopPoint <- all(nextNodes <= tips)
    }
    tipClade <- storeNodes[which(storeNodes <= Ntip(phylo))]
    tippo <- phylo$tip.label[tipClade]
    return(tippo)
}
