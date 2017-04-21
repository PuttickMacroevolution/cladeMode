#' randomNode
#'
#' randomNode
#' @param phylo input phylogeny
#' @param size minimum n tips of descending node (as given by a proportion of total tips)
#' @param overSize return node ancestral to n tips as specified by 'size' or not (default = TRUE)
#' @keywords 
#' @return
#' @export
#' @examples
#' nodeDes()

randomNode <- function(phylo, size=0.25, overSize=TRUE) {
	nodes_phylo <- getNodes(phylo)
	nodes_phylo_0.25 <- ceiling(Ntip(phylo) * size)
	if(overSize == TRUE) randomNode <- nodes_phylo[which(sapply(nodes_phylo, function(x) length(nodeDes(phylo, x))) >= nodes_phylo_0.25)]
	if(overSize != TRUE) randomNode <- nodes_phylo[which(sapply(nodes_phylo, function(x) length(nodeDes(phylo, x))) < nodes_phylo_0.25)]
	return(randomNode)
}
