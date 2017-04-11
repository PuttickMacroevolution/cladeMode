#' convert treeSim to APE
#'
#' convert treeSim to APE
#' @param treeSimPhy phylogeny from treeSim
#' @keywords 
#' @return trApe tree in APE format
#' @export
#' @examples
#' phyIC()

treeSimToApe <- function(treeSimPhy) {
  trApe <- treeSimPhy
  trApe$edge[which(trApe$edge[,2] <= Ntip(trApe)), 2] <- 1:Ntip(trApe)
  original <- c(trApe$edge[1,1], trApe$edge[which(trApe$edge[,2] > Ntip(trApe)), 2]) 
  new <- (Ntip(trApe) + 1) : (Ntip(trApe) +  Nnode(trApe) )
  trApe$edge[which(trApe$edge[,2] > Ntip(trApe)),2] <- new[-1]
  edgeOne <- treeSimPhy$edge[,1]
  
  for(i in 1:length(original)) {
    mover <- which(edgeOne == original[i])
    trApe$edge[mover,1] <- new[i]	
  }
  trApe <- as.phylo(trApe)
  trApe
}