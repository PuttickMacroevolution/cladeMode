#' Simulate nested modes of evolution
#'
#' Simulate various nested modes of continous trait evolution applied to the whole phylogeny (BM, OU, EB) or nested within subclades of the phylogeny ("nested EB", "nested EB rate", "rate Shift", "nested OU")
#' @param n number of datasets to simulate (default = 1)
#' @param phy dated tree in APE format. Although non-ultrametric phylogenies can be used the method has only been tested on ultrametric trees
#' @param mode of evolution for the model. One of "BM", "EB", "OU", "nestedEB", "nestedEBRate", "rateShift", "nestedOU"). "BM" fits a simple model of Brownian Motion; "EB" fits the early burst model (Harmon et al. 2010); "OU" fits the Orstein-Uhlenbeck model (Hansen). The 'nested' models allow shifts to occur in subclades of the phylogeny with an ancestral BM model. 'nested EB' allows for a EB model to a subclade; 'nested EB rate' fits a model a similar model but with the addition of a scalar that allows for an increase on the ancestral BM rate; "nested OU" fits a nested model of the OU process which optimises the alpha parameter; and the "rateShift" model allows for a scalar-based increase or decrease rate in a nested subclade (O'Meara et al. 2006; Thomas et al. 2006)
#' @param node number of node at which to estimate nested evolutionary mode shift
#' @param beta the value of the Brownian motion rate (default = 1) of the ancestral process
#' @param a the value of the early burst parameter only applicable to the 'EB', 'nestedEB', and 'nestedEBRate' models. This value must be smaller than -1e-6, and a lower value signifies a larger exponential decrease through time. If NULL the value is the midpoint of the expected early burst values. 
#' @param shiftRate the value of the shift scalar parameter (default = 1) that signifies an increase (> 1) or decrease (< 1) in the subclade rate of evolution. Only applicable to the 'rateShift' model
#' @param ebBeta the value of the scaler for the nested early burst model (default = 2) that causes an increase in the nested clade alongside the early burst process. If ebBeta=1 this process is identical to the 'nestedEB' model
#' @param alpha the value of the attraction parameter alpha only applicable for the 'OU' and 'nestedOU' models
#' @keywords 
#' @return y a matrix of simulated values
#' @export
#' @examples
#' library(geiger)
#' set.seed(30)
#' # simulate tree with 20 species
#' simTree <- rcoal(20)
#' #Simulate BM data
#' simData <- cladeModeSim(n=1, phy=simTree, mode="BM")
#' #Simulate nested EB with a scalar increase of 5 and a EB rate 'a' of -2
#' simData <- cladeModeSim(n=1, phy=simTree, mode="nestedEB", ebBeta=5, a=-2)

cladeModeSim <- function(n = 1, phy, node=NULL, mode = c("BM", "EB", "nestedEB", "nestedEBRate", "rateShift", "nestedOU", "OU"), beta = 1, a = NULL, alpha = 0.05, shiftRate = 0.5, ebBeta = 2) {
    
    model <- mode
    beta <- log(beta)
    shiftRate <- log(shiftRate)
    ebBeta <- log(ebBeta)
    
    if (model == "BM") {
	        phyMat <- vcv.phylo(phy) * exp(beta)
    	    attr(phyMat, "class") <- "matrix"
    	    y <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
    	    rownames(y) <- phy$tip.label
    	    return(y)
    }
    if (model == "nestedEB") {
    	yy <- c()
        for(x in 1:length(node)) {
	        allTimes <- BranchingTimesFossil(phy)
	        relations_num <- monoClade(phy, node[x])
	        originTimeFocal <- allTimes[which(names(allTimes) == 
            node[x])]
	        names(allTimes)[((Nnode(phy) + 1):(Nnode(phy) + Ntip(phy)))] <- 1:Ntip(phy)
	        times <- c()
	        for (u in 1:length(allTimes)) times[which(phy$edge[, 
            1] == names(allTimes)[u])] <- allTimes[u]
	        branches <- match(relations_num, phy$edge[, 2])
   		     if (is.null(a) == T) 
            a <- (log(1e-04)/originTimeFocal)/2
        	release.mat <- split.vcv.node(phy, node[x])
        	bm.mat <- release.mat[[1]]
        	eb.mat <- ebSubClade(phy, a = a, node = node[x], branches, 
            times, originTimeFocal)
        	vv <- exp(beta) * (bm.mat + eb.mat)
        	y <- as.matrix(t(rmvnorm(1, sigma = vv)))
        	yy <- cbind(yy, y)
        	}
        	
        y <- yy
        rownames(y) <- phy$tip.label
        return(y)
    }
    if (model == "nestedEBRate") {
    	yy <- c()
        for(x in 1:length(node)) {
        	allTimes <- BranchingTimesFossil(phy)
        	relations_num <- monoClade(phy, node[x])
        	originTimeFocal <- allTimes[which(names(allTimes) == 
            node[x])]
        	names(allTimes)[((Nnode(phy) + 1):(Nnode(phy) + Ntip(phy)))] <- 1:Ntip(phy)
        	times <- c()
        	for (u in 1:length(allTimes)) times[which(phy$edge[, 
            1] == names(allTimes)[u])] <- allTimes[u]
        	branches <- match(relations_num, phy$edge[, 2])
        	if (is.null(a) == T) 
            a <- (log(1e-04)/originTimeFocal)/2
        	release.mat <- split.vcv.node(phy, node[x])
        	bm.mat <- release.mat[[1]]
        	eb.mat <- ebSubClade(phy, a = a, node[x], branches, times, 
            originTimeFocal)
        	vv <- exp(beta) * (bm.mat + (exp(ebBeta) * eb.mat))
        	y <- as.matrix(t(rmvnorm(1, sigma = vv)))
        	yy <- cbind(yy, y)
        	}
        y <- yy
       	rownames(y) <- phy$tip.label
       	return(y)
    }
    if (model == "rateShift") {
    	yy <- c()
        for(x in 1:length(node)) {
        	release.mat <- split.vcv.node(phy, node[x])
        	bm.mat <- release.mat[[1]]
        	shift.mat <- release.mat[[2]]
        	vv <- exp(beta) * (bm.mat + (exp(shiftRate) * shift.mat))
        	y <- as.matrix(t(rmvnorm(1, sigma = vv)))
        	yy <- cbind(yy, y)
        	}
        y <- yy
        rownames(y) <- phy$tip.label
        return(y)
    }
    if (model == "nestedOU") {
    	yy <- c()
        for(x in 1:length(node)) {
    	    release.mat <- split.vcv.node(phy, node[x])
        	bm.mat <- release.mat[[1]]
        	shift.mat <- release.mat[[2]]
        	vv <- exp(beta) * (bm.mat + ouMatrix(shift.mat, exp(alpha)))
        	y <- as.matrix(t(rmvnorm(1, sigma = vv)))
        	yy <- cbind(yy, y)
        	}
        y <- yy
        rownames(y) <- phy$tip.label
        return(y)
    }
    if (model == "EB") {
        node = Ntip(phy) + 1
        allTimes <- BranchingTimesFossil(phy)
        relations_num <- monoClade(phy, node)
        originTimeFocal <- allTimes[which(names(allTimes) == 
            node)]
        names(allTimes)[((Nnode(phy) + 1):(Nnode(phy) + Ntip(phy)))] <- 1:Ntip(phy)
        times <- c()
        for (u in 1:length(allTimes)) times[which(phy$edge[, 
            1] == names(allTimes)[u])] <- allTimes[u]
        if (is.null(a) == T) 
            a <- (log(1e-04)/originTimeFocal)/2
        branches <- match(relations_num, phy$edge[, 2])
        eb.mat <- ebSubClade(phy, a = a, node = node, branches = branches, 
            times = times, originTime = originTimeFocal)
        vv <- exp(beta) * (eb.mat)
        y <- as.matrix(t(rmvnorm(n, sigma = vv)))
        rownames(y) <- phy$tip.label
        return(y)
    }
    if (model == "OU") {
        vcvMat <- vcv(phy)
        vv <- exp(beta) * ouMatrix(vcvMat, exp(alpha))
        y <- as.matrix(t(rmvnorm(n, sigma = vv)))
        rownames(y) <- phy$tip.label
        return(y)
    }
}