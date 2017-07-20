#' reScalePhylo function
#'
#' reScalePhylo (internal function)
#' @param phy
#' @param y
#' @param param
#' @param model
#' @keywords 
#' @return list phylogeny and matrix of parameters
#' @export
#' @examples
#' reScalePhylo()

reScalePhylo <- function(phy, y, param, model=c("BM", "EB", "nestedEB", "nestedEBRate", "rateShift", "rateSlow", "nestedOU")) {
				
	if(model == "BM") {
	 	y <- y[match(phy$tip.label, names(y))]
		sig <- log(param[3])
		phy$edge.length <- phy$edge.length * exp(sig)
		return(list(phy, y))
	}
	
	if(model == "nestedEB") {
		sig <- log(param[3])
		a <- param[4]
		node <- param[5]
		y <- y[match(phy$tip.label, names(y))]
 		allTimes <- branching.times(phy)
		relations_num <- monoClade(phy, node)
  		originTimeFocal <- allTimes[which(names(allTimes) == node)]
  		times <- c()
  		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  		branches <- match(relations_num , phy$edge[,2])
  		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
  		phy <- subClade(phy, a=a, node=node, branches=branches, times=times, originTime=originTimeFocal)
		phy$edge.length <- phy$edge.length * exp(sig)
		return(list(phy, y))
	}
			
	if(model == "nestedEBRate") {
		sig <- log(param[3])
		a <- param[4]
		node <- param[6]
		ebRate <- log(param[5])
		y <- y[match(phy$tip.label, names(y))]
 		allTimes <- branching.times(phy)
  		relations_num <- monoClade(phy, node)
  		originTimeFocal <- allTimes[which(names(allTimes) == node)]
  		times <- c()
  		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  		branches <- match(relations_num , phy$edge[,2])
  		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
  		phy <- subClade(phy, a=a, node=node, branches=branches, times=times, originTime=originTimeFocal, rate=exp(ebRate))
		phy$edge.length <- phy$edge.length * exp(sig)
		return(list(phy, y))

		}
  		  		
  	if(model == "rateShift") {
		sig <- log(param[3])
		node <- param[5]
		rate <- log(param[4])
		y <- y[match(phy$tip.label, names(y))]
 		allTimes <- branching.times(phy)
  		relations_num <- monoClade(phy, node)
  		originTimeFocal <- allTimes[which(names(allTimes) == node)]
  		times <- c()
  		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  		branches <- match(relations_num , phy$edge[,2])
  		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
		phy <- rateClade(phy, branches=branches, rate=exp(rate))
		phy$edge.length <- phy$edge.length * exp(sig)
		return(list(phy, y))
  		}
  		  		
 	if(model == "nestedOU") {
		sig <- log(param[3])
		node <- param[5]
		alpha <- log(param[4])
		y <- y[match(phy$tip.label, names(y))]
 		allTimes <- branching.times(phy)
		originTimeFocal <- allTimes[which(names(allTimes) == node)]
		relations_num <- monoClade(phy, node)
  		times <- c()
  		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  		branches <- match(relations_num , phy$edge[,2])
  		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]		
  		phy <- OUtree(phy, branches=branches, times=times, alpha=exp(alpha), originTime=originTimeFocal)
		phy$edge.length <- phy$edge.length * exp(sig)
		return(list(phy, y))
  		}
  			
  		
	if(model == "EB") {
		sig <- log(param[3])
		a <- param[4]
		y <- y[match(phy$tip.label, names(y))]
 		allTimes <- branching.times(phy)
		node=Ntip(phy) + 1
		relations_num <- monoClade(phy, node)
		originTimeFocal <- allTimes[which(names(allTimes) == node)]
		times <- c()
		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
		branches <- match(relations_num, phy$edge[,2])
		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
 		phy <- subClade(phy, a=a, node=node, branches=branches, times=times, originTime=originTimeFocal)
		phy$edge.length <- phy$edge.length * exp(sig)
		return(list(phy, y))
		}
  	
		
	if(model == "OU") {
		sig <- log(param[3])
		alpha <- log(param[4])
		y <- y[match(phy$tip.label, names(y))]
 		allTimes <- branching.times(phy)
		node=Ntip(phy) + 1
		relations_num <- monoClade(phy, node)
		originTimeFocal <- allTimes[which(names(allTimes) == node)]
		times <- c()
		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
		branches <- match(relations_num, phy$edge[,2])
		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
		phy <- OUtree(phy, alpha=exp(alpha), branches=branches, times=times, originTime=originTimeFocal)
		phy$edge.length <- phy$edge.length * exp(sig)
		return(list(phy, y))
 		}

}
