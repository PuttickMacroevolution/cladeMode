#' Analyse nested modes of evolution
#'
#' Analyses various nested modes of continous trait evolution applied to the whole phylogeny (BM, OU, EB) or nested within subclades of the phylogeny ("nested EB", "nested EB rate", "rate Shift", "nested OU")
#' @param phy dated tree in APE format. Although non-ultrametric phylogenies can be used the method has only been tested on ultrametric trees
#' @param y named vector of trait data
#' @param node number of node at which to estimate nested evolutionary mode. If FALSE (default) the model tests model fit at all nodes using an iterative procedure. The nodes at which the model is applied it determined by the minSize argument
#' @param minSize the minimum size of node at which to select a nested rate. The number represents the minimum size of the clade for which all models will be nested. 
#' @param mode of evolution for the model. One of "BM", "EB", "OU", "nestedEB", "nestedEBRate", "rateShift", "nestedOU"). "BM" fits a simple model of Brownian Motion; "EB" fits the early burst model (Harmon et al. 2010); "OU" fits the Orstein-Uhlenbeck model (Hansen). The 'nested' models allow shifts to occur in subclades of the phylogeny with an ancestral BM model. 'nested EB' allows for a EB model to a subclade; 'nested EB rate' fits a model a similar model but with the addition of a scalar that allows for an increase on the ancestral BM rate; "nested OU" fits a nested model of the OU process which optimises the alpha parameter; and the "rateShift" model allows for a scalar-based increase or decrease rate in a nested subclade (O'Meara et al. 2006; Thomas et al. 2006)
#' @param ncores number of seperate cores to simultaneously estimate rate shifts - defaults to 1
#' @param contrasts logical (default = FALSE) as to whether to estimate parameters with the default likelihood search or independent contrasts
#' @param returnPhy logical (default = FALSE) returns re-scaled phylogeny according the best-fitting model
#' @param measureError a named vector of measurement error around the trait data y. If NULL (default) the model assumes no measurement error. Warning: this is currently not fully tested and only applies to all models (except BM) when 'contrasts=FALSE'
#' @keywords 
#' @return results a list containing the model parameters, log likelihood, aicc, optionally the rescaled tree
#' @export
#' @examples
#' library(geiger)
#' set.seed(30)
#' # simulate tree with 20 species
#' simTree <- rcoal(20)
#' simData <- cladeModeSim(n=1, phy=simTree, mode="BM")
#' ##### comparison between geiger and cladeMode parameter estimates using early burst model
#' #cladeMode
#' unlist(cladeMode(phy=simTree, y=simData[,1], mode="EB", cont=T))
#' #cladeMode paramater estimates 
#' #         lnl   root.state         beta            a          aic         aicc            k 
#' #-13.12825138  -0.38104635   1.14123829  -0.06945831  32.25650276  33.85650276   3.00000000 
# '#fitContinuous
#' unlist(fitContinuous(simTree, simData[,1], model="EB")[[4]])
#' # a          sigsq       z0           lnL        method	k	aic			aicc
#' # -0.06945   1.14123  -0.3810463  -13.1282513    subplex	3	32.2565		33.7565
#' # The parameters are identical but there is one discrepency - cladeMode uses n in the aicc calculation as the number of taxa (i.e, Butler and King 2004; O'Meara et al. 2006), but Geiger uses n + 1 . This should not affect results as only cladeMode is used in all analyses and thus comparison of AICc values are consistent, but care should be taken if compared against Geiger models
#' # cladeMode AICc
#' n=20
#' k=3
#' 2 * k * (n - 1)/(n - k - 2) - 2 * -13.12825138
#' # Geiger AICc
#' n=21
#' k=3
#' 2 * k * (n - 1)/(n - k - 2) - 2 * -13.12825138
#' ################### Comparison of models            
#' # example use of each model in the function cladeMode          
#' # bm  
#' bmRes <- cladeMode(phy=simTree, y=simData[,1], mode="BM", cont=T)
#' # early burst
#' ebRes <-cladeMode(phy=simTree, y=simData[,1], mode="EB", cont=T)
# nested early burst with concurrent rate increase
#' nestedEBRateRes <-cladeMode(phy=simTree, y=simData[,1], mode="nestedEBRate", cont=T)
#' # nested early burst 
#' nestedEBRes <-cladeMode(phy=simTree, y=simData[,1], mode="nestedEB", cont=T)
#' # nested OU
#' nestedOURes <-cladeMode(phy=simTree, y=simData[,1], mode="nestedOU", cont=T)
#' # nested rate shift
#' nestedRateShiftRes <-cladeMode(phy=simTree, y=simData[,1], mode="rateShift", cont=T)
#' #comparison of models to summarise the best relative fit
#' aicW <- aiccWeight(c(bmRes$aicc, ebRes$aicc, nestedEBRateRes$aicc, nestedEBRes$aicc, nestedOURes$aicc, nestedRateShiftRes$aicc))
#' rownames(aicW) <- c("bm", "eb", "nestedEBRate", "nestedEB", "nestedOU", "nestedRateShift")
#' aicW
#' # nested rate shift has highest, but erroneus support (data were simulated under BM)
#' # Error correction. We will simulate 100 datasets (1000 were used in the ms, 100 used here for brevity) under BM and identify the error rate for the nested rate shift model
#' hunderedSim <- cladeModeSim(n=100,simTree, model="BM")
#' # get aicc of nested rate models
#' nestedRateAICcSim <- sapply(1:100, function(x) cladeMode(phy=simTree, y=hunderedSim[,x], mode="rateShift", cont=T)$aicc)
#' # get aicc of bm models
#' bmAICcSim <- sapply(1:100, function(x) cladeMode(phy=simTree, y=hunderedSim[,x], mode="BM", cont=T)$aicc)
#' #calculate the 95th quantile under which nested rate is wrongly support
#' cutOffAICc <- as.numeric(quantile(nestedRateAICcSim - bmAICcSim, 0.95))
#' # add this number to nested rate aicc score and bm receives true model support
#' aicW <- aiccWeight(c(bmRes$aicc, nestedRateShiftRes$aicc + cutOffAICc))
#' aicW

cladeMode <- function(phy, y, node=F, minSize=5, mode=c("BM", "OU", "EB",  "nestedEB", "nestedEBRate", "rateShift", "nestedOU"), Ncores=1, contrasts=FALSE, returnPhy=FALSE, measureError=NULL) {
	require(parallel)
	require(ape)
	require(mvtnorm)
	
	cont <- contrasts
	nested <- c("nestedEB", "nestedEBRate", "rateShift", "nestedOU")
	if(is.null(measureError)) measureError <- rep(0, Ntip(phy))
	
	mErr <- measureError 
	
	if(sum(node) != 0 && is.na(match(mode, nested))) print("mode applies to whole tree: node argument ignored")
	if(is.null(names(y)) == T) stop("please provide data names")
	if(length(measureError) != Ntip(phy)) stop("measurement error not equal in length to trait data")
	
	
	if(sum(match(mode, nested), na.rm = T) == 0) {
		return(cladeModeNode(phy, y, model=mode, cont=cont, returnPhy=returnPhy, mErr=mErr))
	}
	
	if(is.numeric(node) == F && sum(match(mode, nested), na.rm = T) != 0) {
		allNodesComplete <- getNodes(phy)
		tips <- sapply(1:length(allNodesComplete), function(u) length(which(monoClade(phy, allNodesComplete[u]) <= Ntip(phy))))
		allNodes <- allNodesComplete[which(tips > minSize)] 
		allNest <- mclapply(allNodes, mc.cores=Ncores, function(x) cladeModeNode(phy, y, node=x, model=mode, cont=cont, returnPhy=returnPhy, mErr=mErr))
		converge <- sapply(allNest, function(x) is.numeric(x[[1]]))
		allNest <- allNest[converge]
		aiccScore <- sapply(1:length(allNest), function(x) allNest[[x]]$aicc)
		bestModel <- allNest[which(aiccScore == min(aiccScore))][[1]]
		return(bestModel)
	}

	if(is.numeric(node) == T && sum(match(mode, nested), na.rm = T) != 0) {
		return(cladeModeNode(phy, y, node=node, model=mode, cont=cont, returnPhy=returnPhy, mErr=mErr))
	}
}