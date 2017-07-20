#' rescalePIC function
#'
#' rescalePIC a function to produce a re-scaled phylogeny using IC
#' @param phy dated tree in ape format
#' @param y named vector of trait values
#' @keywords 
#' @return phy a rescaled phylogeny
#' @export
#' @examples
#' rescalePIC()

rescalePIC <- function(phy, y) {
	bL <- phy$edge.length
	tips <- which(phy$edge[,2] <= Ntip(phy))
	tipsMat <- phy$edge[tips,]
	n <- Ntip(phy)
	rescaled <- bL
	rescaled[tips] <- bL[tips]
	internal <- which(phy$edge[,2] > Ntip(phy))
	allInt <- rev(c(Ntip(phy)+1, phy$edge[internal, 2]))

	for(u in 1:length(allInt)) {
		nodals <- which(allInt[u] == phy$edge[,1])
		resc <- rescaled[nodals]
		if(sum(which(resc == 0)) != 0) resc <- resc[-which(resc == 0)]
		addOn <- prod(resc) / sum(resc)
		toL <- which(allInt[u] == phy$edge[,2])
		rescaled[toL] <- rescaled[toL] + addOn
	}
	
	phy$edge.length <- rescaled
	return(phy)
	
}

