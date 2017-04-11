#' Simulate nested modes of evolution
#'
#' Simulate various nested modes of trait evolution
#' @param n number of datasets to simulate (default = 1)
#' @param phy tree in APE formate
#' @param y named vector of trait data
#' @param node number of node at which to estimate nested evolutionary mode
#' @param mode of evolution for the model. One of "BM", "EB", "nestedEB", "nestedEBRate", "rateShift", "nestedOu", "rateSlow")
#' @param beta the value of the Brownian motion rate (default = 1) 
#' @param a the value of the early burst parameter. If NULL the value is the midpoint of the expected early burst values
#' @param shiftRate the value of the shiftRate parameter (default = 1) 
#' @param ebBeta the value of the scaler for the nested early burst model (default = 2) 
#' @param mode of evolution for the model. One of "BM", "EB", "nestedEB", "nestedEBRate", "rateShift", "nestedOu", "rateSlow")
#' @keywords 
#' @return y a matrix of simulated values
#' @export
#' @examples
#' cladeMode()


cladeModeSim <- function(n=1, phy, node, model=c("BM", "EB", "nestedEB", "nestedEBRate", "rateShift", "nestedOU", "OU"), beta=1, a=NULL, alpha=0.05, shiftRate=0.5, ebBeta=2) {
	
	beta <- log(beta)
	shiftRate <- log(shiftRate)
	ebBeta <- log(ebBeta)

	if(model == "BM") {
		phyMat <- vcv.phylo(phy) * exp(beta)
		attr(phyMat, "class") <- "matrix"
		y <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
		rownames(y) <- phy$tip.label
		return(y)
	}
	
	if(model == "nestedEB") {
		allTimes <- BranchingTimesFossil(phy)
  		relations_num <- monoClade(phy, node)
  		originTimeFocal <- allTimes[which(names(allTimes) == node)]
  		names(allTimes)[((Nnode(phy)+1):(Nnode(phy)+Ntip(phy)))] <- 1:Ntip(phy)
  		times <- c()
  		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  		branches <- match(relations_num , phy$edge[,2])
   		if(is.null(a) == T) a <- (log(1e-4)/originTimeFocal)	/2	
   		release.mat <- split.vcv.node(phy, node)
  		bm.mat <- release.mat[[1]]	
   		eb.mat <- ebSubClade(phy, a=a, node=node, branches, times, originTimeFocal)
    	vv <- exp(beta) * (bm.mat + eb.mat)
		y <- as.matrix(t(rmvnorm(n,  sigma = vv)))
		rownames(y) <- phy$tip.label
		return(y)
		}
  
	if(model == "nestedEBRate") {
		allTimes <- BranchingTimesFossil(phy)
	  	relations_num <- monoClade(phy, node)
	  	originTimeFocal <- allTimes[which(names(allTimes) == node)]
	  	names(allTimes)[((Nnode(phy)+1):(Nnode(phy)+Ntip(phy)))] <- 1:Ntip(phy)
	  	times <- c()
	  	for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
	  	branches <- match(relations_num , phy$edge[,2])
   		if(is.null(a) == T) a <- (log(1e-4)/originTimeFocal)	/2	
   		release.mat <- split.vcv.node(phy, node)
  		bm.mat <- release.mat[[1]]	 		
   		eb.mat <- ebSubClade(phy, a=a, node, branches, times, originTimeFocal)
   		vv <- exp(beta) * (bm.mat + (exp(ebBeta) * eb.mat))  
	    	y <- as.matrix(t(rmvnorm(n,  sigma = vv)))
   	 	rownames(y) <- phy$tip.label
   	 	return(y)
		}
           
	if(model == "rateShift") {
		release.mat <- split.vcv.node(phy, node)
		bm.mat <- release.mat[[1]]
		shift.mat <- release.mat[[2]]
		vv <- exp(beta) * (bm.mat + (exp(shiftRate) * shift.mat))  
    	y <- as.matrix(t(rmvnorm(n,  sigma = vv)))
    	rownames(y) <- phy$tip.label
    	return(y)
		}
 
 	if(model == "nestedOU") {
		release.mat <- split.vcv.node(phy, node)
		bm.mat <- release.mat[[1]]
		shift.mat <- release.mat[[2]]
		vv <- exp(beta) * (bm.mat + ouMatrix(shift.mat, exp(alpha)))
	    y <- as.matrix(t(rmvnorm(n,  sigma = vv)))
	    rownames(y) <- phy$tip.label
	    return(y)
	 	}
 
	if(model == "EB") {
		node=Ntip(phy) + 1
		allTimes <- BranchingTimesFossil(phy)
		relations_num <- monoClade(phy, node)
		originTimeFocal <- allTimes[which(names(allTimes) == node)]
		names(allTimes)[((Nnode(phy)+1):(Nnode(phy)+Ntip(phy)))] <- 1:Ntip(phy)
		times <- c()
		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
		if(is.null(a) == T) a <- (log(1e-4)/originTimeFocal)	/2			
		branches <- match(relations_num, phy$edge[,2])
	   	eb.mat <- ebSubClade(phy, a=a, node=node, branches=branches, times=times, originTime=originTimeFocal)
    		vv <- exp(beta) * (eb.mat)  
    		y <- as.matrix(t(rmvnorm(n,  sigma = vv)))
    		rownames(y) <- phy$tip.label
    		return(y)
		}   
	
	if(model == "OU") {
		vcvMat <- vcv(phy)
	   	vv <- exp(beta) * ouMatrix(vcvMat, exp(alpha))
		y <- as.matrix(t(rmvnorm(n,  sigma = vv)))
		rownames(y) <- phy$tip.label
		return(y)
		}
}