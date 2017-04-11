#' findNodes function
#'
#' findNodes function
#' @param phylo
#' @param minSize (default = 0.1)
#' @keywords 
#' @return node
#' @export
#' @examples
#' findNodes()

findNodes <- function(phylo, minSize=0.1) {
    allNodes <- getNodes(phylo) 
    tips <- sapply(1:length(allNodes), function(u) length(which(monoClade(phylo, allNodes[u]) <= Ntip(phylo))))
    allNodes <- allNodes[which(tips > ceiling(Ntip(phylo) * minSize))]
     if(length(allNodes) > 1) {
      node <- sample(c(allNodes), 1)
    } else
    {
      node <- allNodes
    }
   return(node)
 }
