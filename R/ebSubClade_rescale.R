#' Rescale phylogeny with early burst values
#'
#' Internal function to rescale phylogeny with early burst values
#' @param phyIn input phylogeny
#' @param a Input
#' @param node_in
#' @keywords 
#' @return rescaled vcv matrix
#' @export
#' @examples
#' cladeMode()


ebSubClade_rescale <- function (phyIn, aInput, node_in) {
  
  allTimes <- BranchingTimesFossil(phyIn)
  relations_num <- monoClade(phy=phyIn, node_in)
  originTime <- allTimes[which(names(allTimes) == node_in)]
  names(allTimes)[((Nnode(phyIn)+1):(Nnode(phyIn)+Ntip(phyIn)))] <- 1:Ntip(phyIn)
  times <- c()
  for(u in 1:length(allTimes)) times[which(phyIn$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  branches <- match(relations_num , phyIn$edge[,2])
  phyEB <- phyIn
  
  for (i in branches) {
    bl <- phyIn$edge.length[i]
    age = times[i]
    t1 = (originTime - age) 
    t2 = (t1 + bl) 
    phyEB$edge.length[i] = (exp(aInput*t2)-exp(aInput*t1))/(aInput)
  }
  return(phyEB)
}