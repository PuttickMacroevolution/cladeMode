#' Analyse nested modes of evolution
#'
#' Analyses various nested modes of trait evolution
#' @param phy tree in APE formate
#' @param y named vector of trait data
#' @param node number of node at which to estimate nested evolutionary mode. If FALSE (default) the model searches across all nodes
#' @param minSize the minimum size of node at which to select a nested rate
#' @param mode of evolution for the model. One of "BM", "EB", "nestedEB", "nestedEBRate", "rateShift", "nestedOu", "rateSlow")
#' @param ncores number of seperate cores to simultaneously estimate rate shifts - defaults to 1
#' @param contrasts logical (default = FALSE) as to whether to estimate parameters with likelihood search or independent contrasts
#' @param returnPhy logical (default = FALSE) return re-scaled phylogeny
#' @keywords 
#' @return results a list containing the model parameters, log likelihood, aicc, optionally the rescaled tree
#' @export
#' @examples
#' cladeMode()

cladeMode <- function(phy, y, node=F, minSize=5, mode=c("BM", "EB",  "nestedEB", "nestedEBRate", "rateShift", "nestedOU", "rateSlow"), Ncores=1, contrasts=F, returnPhy=FALSE) {
	require(parallel)
	require(ape)
	require(mvtnorm)
	
	cont <- contrasts
	nested <- c("nestedEB", "nestedEBRate", "rateShift", "nestedOU", "rateSlow")

	if(sum(node) != 0 && is.na(match(mode, nested))) print("mode applies to whole tree: node argument ignored")
	if(is.null(names(y)) == T) stop("please provide data names")
	
	if(sum(match(mode, nested), na.rm = T) == 0) {
		return(cladeModeNode(phy, y, model=mode, cont=cont, returnPhy=returnPhy))
	}
	
	if(is.numeric(node) == F && sum(match(mode, nested), na.rm = T) != 0) {
		allNodesComplete <- getNodes(phy)
		tips <- sapply(1:length(allNodesComplete), function(u) length(which(monoClade(phy, allNodesComplete[u]) <= Ntip(phy))))
		allNodes <- allNodesComplete[which(tips > minSize)] 
		allNest <- mclapply(allNodes, mc.cores=Ncores, function(x) cladeModeNode(phy, y, node=x, model=mode, cont=cont, returnPhy=returnPhy))
		converge <- sapply(allNest, function(x) is.numeric(x[[1]]))
		allNest <- allNest[converge]
		aiccScore <- sapply(1:length(allNest), function(x) allNest[[x]]$aicc)
		bestModel <- allNest[which(aiccScore == min(aiccScore))][[1]]
		return(bestModel)
	}

	if(is.numeric(node) == T && sum(match(mode, nested), na.rm = T) != 0) {
		return(cladeModeNode(phy, y, node=node, model=mode, cont=cont, returnPhy=returnPhy))
	}
}