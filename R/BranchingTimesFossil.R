#' BranchingTimesFossil function from Graham Slater
#'
#' BranchingTimesFossil function from Graham Slater
#' @param phy
#' @keywords 
#' @return list phylogeny and matrix of parameters
#' @export
#' @examples
#' phylogMean()

BranchingTimesFossil <- function (phy, tol = .Machine$double.eps^0.5) 
{
    if (class(phy) != "phylo") 
        stop("object \"phy\" is not of class \"phylo\"")
    phy2 <- phy
    phy <- new2old.phylo(phy)
    tmp <- as.numeric(phy$edge)
    nb.tip <- max(tmp)
    nb.node <- -min(tmp)
    xx <- as.numeric(rep(NA, nb.tip + nb.node))
    names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
    xx["-1"] <- 0
    for (i in 2:length(xx)) {
        nod <- names(xx[i])
        ind <- which(as.numeric(phy$edge[, 2]) == nod)
        base <- phy$edge[ind, 1]
        xx[i] <- xx[base] + phy$edge.length[ind]
    }
    bt <- abs(xx - max(xx));
	
	for(i in 1:length(bt)) {
		
		if(bt[i]<.Machine$double.eps^0.5) bt[i] <- 0; 	}
	
	names(bt) <- c(seq(nb.tip+1, nb.tip+nb.node), phy$tip.label)
	
	
	return(bt);
}
