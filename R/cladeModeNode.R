#' Analyse nested modes of evolution - internal function
#'
#' Internal function to analyses various nested modes of trait evolution
#' @param 
#' @keywords 
#' @return
#' @export
#' @examples
#' 

cladeModeNode <- function(phy, y, node, model=c("BM", "EB", "nestedEB", "nestedEBRate", "rateShift", "nestedOU"), cont, returnPhy, mErr) {
	
	IC <- cont
 	allTimes <- branching.times(phy)
 	y <- y[match(phy$tip.label, names(y))]
    n <- Ntip(phy)
    beta.start <- var(y)/max(allTimes)

    startBM = log(beta.start)
    lowerBM = log(1e-8)
    upperBM = log(20)
	k <- 2
	picMark <- phyIC(phy, y)
	bmLik <- picMark$logLikelihood
			
	if(model == "BM") {
		results <- list(lnl = picMark$logLikelihood, root.state = as.numeric(picMark$phyloMean), beta = picMark$sigmaSq)
		results$aic <- 2 * k - 2 * picMark$logLikelihood
  		results$aicc <- 2 * k * (n - 1)/(n - k - 2) - 2 * picMark$logLikelihood
   		results$aicc <- aicc(picMark$logLikelihood, k, Ntip(phy))	
  		results$k <- k
  		if(returnPhy == TRUE) results$phy <- phy
  		return(results)
	}
	
	if(model == "nestedEB") {
		relations_num <- monoClade(phy, node)
  		originTimeFocal <- allTimes[which(names(allTimes) == node)]
  		times <- c()
  		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  		branches <- match(relations_num , phy$edge[,2])
  		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
  		k <- 3	
  		
  		if(IC == F) {
  			start <- c(startBM, -0.01)
  			upper <- c(upperBM, -0.000001)
  			lower <- c(lowerBM, ACDC.prior(phy)[[1]])
  			bm.mat <- split.vcv.node(phy, node)[[1]]
			foo <- function(x) {
    			eb.mat <- ebSubClade(phy, a=x[2], node=node, branches, times, originTimeFocal)
    			vv <- exp(x[1]) * (bm.mat + eb.mat)
    			diag(vv) <- diag(vv) + mErr ^ 2
    			mu <- phylogMean(vv, y)  
     			mu <- rep(mu, n)
      			return(-dmvnorm(y, mu, vv, log = T))
   				}
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
			ml.vcv <- exp(o$par[1]) * (bm.mat + ebSubClade(phy, a=o$par[2], node=node, branches=branches, times=times, originTime=originTimeFocal))
 			root.state <- phylogMean(ml.vcv, y)
 			results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), a = o$par[2], focalnode = node);
 			if(returnPhy == TRUE) results$phy <- phy

			} else {
				start <- -0.01
  				upper <- -0.000001
  				lower <- log(10^-5)/originTimeFocal
				foo <- function(x) {
					eb.mat <- subClade(phy, a=x, node=node, branches=branches, times=times, originTime=originTimeFocal)
					return(-phyIC(eb.mat, y)$logLikelihood)
					}
				o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
				mlPhy <- subClade(phy, a=o$par, node=node, branches=branches, times=times, originTime=originTimeFocal)
				xx <- phyIC(mlPhy, y)
				results <- list(lnl = xx$logLikelihood, root.state = as.numeric(xx$phyloMean), beta = xx$sigmaSq, a = o$par, focalnode=node);
				results$aic <- 2 * k - 2 * results$lnl
				results$aicc <- 2 * k * (n - 1)/(n - k - 2) - 2 * results$lnl
				if(returnPhy == TRUE) results$phy <- xx$rescaledPhylo
			}
			
			results$aic <- 2 * k - 2 * results$lnl
			results$aicc <- aicc(results$lnl, k, Ntip(phy))
			results$k <- k
			return(results)
    }

	if(model == "nestedEBRate") {

  		relations_num <- monoClade(phy, node)
  		originTimeFocal <- allTimes[which(names(allTimes) == node)]
  		times <- c()
  		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  		branches <- match(relations_num , phy$edge[,2])
  		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
  		k <- 4
  		if(IC == F) {
			start <- c(startBM, log(2),  -0.01)
			upper <- c(upperBM, upperBM, -0.000001)
			lower <- c(lowerBM, log(1),  ACDC.prior(phy)[[1]])
			bm.mat <- split.vcv.node(phy, node)[[1]]
    		foo <- function(x) {
  		  		eb.mat <- ebSubClade(phy, a=x[3], node, branches, times, originTimeFocal)
    				vv <- exp(x[1]) * (bm.mat + (exp(x[2]) * eb.mat))  
    				diag(vv) <- diag(vv) + mErr ^ 2
    				mu <- phylogMean(vv, y)  
     			mu <- rep(mu, n)
      			return(-dmvnorm(y, mu, vv, log = T))
    			}
 		o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
 		ml.vcv <- exp(o$par[1]) * (bm.mat + (exp(o$par[2]) * ebSubClade(phy, a=o$par[3], node, branches, times, originTimeFocal)))
  		root.state <- phylogMean(ml.vcv, y)
  		results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), a = o$par[3], ebRate = exp(o$par[2]), focalnode = node);
  		} else 
  		{
  			start <- c(log(2),  -0.01)
			upper <- c(upperBM, -0.000001)
			lower <- c(log(1),  log(10^-5)/originTimeFocal)
			foo <- function(x) {
				eb.mat <- subClade(phy, a=x[2], node=node, branches=branches, times=times, originTime=originTimeFocal, rate=exp(x[1]))
				return(-phyIC(eb.mat, y)$logLikelihood)
				}
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
			mlPhy <- subClade(phy, a=o$par[2], node=node, branches=branches, times=times, originTime=originTimeFocal, rate=exp(o$par[1]))
			xx <- phyIC(mlPhy, y)
			results <- list(lnl = xx$logLikelihood, root.state = as.numeric(xx$phyloMean), beta = xx$sigmaSq, a = o$par[2], ebRate = exp(o$par[1]), focalnode = node);
			if(returnPhy == TRUE) results$phy <- xx$rescaledPhylo
		}
  		  		
   		results$aic <- 2 * k - 2 * results$lnl
   		results$aicc <- aicc(results$lnl, k, Ntip(phy))
  		results$k <- k
  		return(results)
	}

	if(model == "rateShift") {
	
		k <- 3
  		if(IC == F) {
			start <- c(startBM, startBM)
			upper <- c(upperBM, upperBM)
			lower <- c(lowerBM, lowerBM)
			release.mat <- split.vcv.node(phy, node)
			bm.mat <- release.mat[[1]]
			shift.mat <- release.mat[[2]]
 			foo <- function(x) {
				vv <- exp(x[1]) * (bm.mat + (exp(x[2]) * shift.mat))  
    			diag(vv) <- diag(vv) + mErr ^ 2
    			mu <- phylogMean(vv, y)  
    		  	mu<-rep(mu, n)
      			return(-dmvnorm(y, mu, vv, log = T))
    			}
  
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
  			ml.vcv <- exp(o$par[2]) * (bm.mat + (exp(o$par[1]) * shift.mat))
  			root.state <- phylogMean(ml.vcv, y)
  			results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), shiftRate = exp(o$par[2]), focalnode = node);
  		} else 
  		{
  			relations_num <- monoClade(phy, node)
  			originTimeFocal <- allTimes[which(names(allTimes) == node)]
  			times <- c()
  			for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  			branches <- match(relations_num , phy$edge[,2])
  			if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
			start <- startBM
			upper <- upperBM
			lower <- lowerBM
			foo <- function(x) {
				eb.mat <- rateClade(phy, branches=branches, rate=exp(x))
				return(-phyIC(eb.mat, y)$logLikelihood)
				}
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
			mlPhy <- rateClade(phy, branches=branches, rate=exp(o$par))
			xx <- phyIC(mlPhy, y)
			results <- list(lnl = xx$logLikelihood, root.state = as.numeric(xx$phyloMean), beta = xx$sigmaSq, shiftRate = exp(o$par), focalnode = node);
			if(returnPhy == TRUE) results$phy <- xx$rescaledPhylo

  		}
  		  		
  		results$aic <- 2 * k - 2 * results$lnl
		results$aicc <- aicc(results$lnl, k, Ntip(phy))
  		results$k <- k
  		return(results)
  	}
  	  	
 	if(model == "nestedOU") {

		k <- 3
		if(IC == F){
			start <- c(startBM, log(0.05))
			upper <- c(upperBM, 1)
			lower <- c(lowerBM, log(1e-8))
			release.mat <- split.vcv.node(phy, node)
			n <- Ntip(phy)
			bm.mat <- release.mat[[1]]
			shift.mat <- release.mat[[2]]
  
			foo <- function(x) {
    			vv <- exp(x[1]) * (bm.mat + ouMatrix(shift.mat, exp(x[2])))
				diag(vv) <- diag(vv) + mErr ^ 2
				mu <- phylogMean(vv, y)  
      			mu<-rep(mu, n)
      			return(-dmvnorm(y, mu, vv, log = T))
    			}
  
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
	  		ml.vcv <- exp(o$par[1]) * (bm.mat + (ouMatrix(shift.mat, exp(o$par[2]))))
	  		root.state <- phylogMean(ml.vcv, y)
	  		results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]), focalnode = node);
  		} else 
  		{
  			relations_num <- monoClade(phy, node)
  			originTimeFocal <- allTimes[which(names(allTimes) == node)]
  			times <- c()
  			for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
  			branches <- match(relations_num , phy$edge[,2])
  			if(any(is.na(branches))) branches <- branches[complete.cases(branches)]		
  			start <- log(0.05)
			upper <- log(exp(1))
			lower <- log(1e-8)
  			foo <- function(x) {
  				tr <- OUtree(phy, exp(x), branches=branches, times, originTimeFocal)
				return(-phyIC(tr, y)$logLikelihood)
				}
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
			mlPhy <- OUtree(phy, branches=branches, times=times, alpha=exp(o$par), originTime=originTimeFocal)
			xx <- phyIC(mlPhy, y)
			results <- list(lnl = xx$logLikelihood, root.state = as.numeric(xx$phyloMean), beta = xx$sigmaSq, alpha = exp(o$par), focalnode = node);
			if(returnPhy == TRUE) results$phy <- xx$rescaledPhylo
  		}
  			
  			results$aic <- 2 * k - 2 * results$lnl
			results$aicc <- aicc(results$lnl, k, Ntip(phy))
  			results$k <- k
  			return(results)  
  		}

	if(model == "EB") {
		
		
		node=Ntip(phy) + 1
		relations_num <- monoClade(phy, node)
		originTimeFocal <- allTimes[which(names(allTimes) == node)]
		times <- c()
		for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
		branches <- match(relations_num, phy$edge[,2])
		if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
 		k <- 3
 		if(IC == F)	{
			start <- c(startBM, -0.01)
			upper <- c(upperBM, -0.000001)
			lower <- c(lowerBM, ACDC.prior(phy)[[1]])
	  		foo <- function(x) {
    			eb.mat <- ebSubClade(phy, a=x[2], node=node, branches=branches, times=times, originTime=originTimeFocal)
    			vv <- exp(x[1]) * (eb.mat)  
    			diag(vv) <- diag(vv) + mErr ^ 2
    			mu <- phylogMean(vv, y)  
      			mu<-rep(mu, n)
      			return(-dmvnorm(y, mu, vv, log = T))
    		}
   			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
			ml.vcv <- exp(o$par[1]) * ebSubClade(phy, a=o$par[2], node=node, branches=branches, times=times, originTime=originTimeFocal)
  			root.state <- phylogMean(ml.vcv, y)
  			results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), a = o$par[2]);
  		} else
  		{
  			start <- -0.01
			upper <- -0.000001
			lower <- log(10^-5)/originTimeFocal
					
			foo <- function(x) {
				eb.mat <- subClade(phy, a=x, node=node, branches=branches, times=times, originTime=originTimeFocal)
				return(-phyIC(eb.mat, y)$logLikelihood)
				}
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
			mlPhy <- subClade(phy, a=o$par, node=node, branches=branches, times=times, originTime=originTimeFocal)
			xx <- phyIC(mlPhy, y)
			results <- list(lnl = xx$logLikelihood, root.state = as.numeric(xx$phyloMean), beta = xx$sigmaSq, a = o$par);
			if(returnPhy == TRUE) results$phy <- xx$rescaledPhylo
  		}
  		  		
  		results$aic <- 2 * k - 2 * results$lnl
		results$aicc <- aicc(results$lnl, k, Ntip(phy))
  		results$k <- k
  		return(results) 
	}
		
	if(model == "OU") {

		k <- 3
		if(IC == F){
			 k <- 3
        start <- c(startBM, log(0.05))
        upper <- c(upperBM, log(exp(1)))
        lower <- c(lowerBM, log(1e-08))
        vcvMat <- vcv(phy)
        foo <- function(x) {
            vv <- exp(x[1]) * ouMatrix(vcvMat, exp(x[2]))
           	diag(vv) <- diag(vv) + mErr ^ 2
            mu <- phylogMean(vv, y)
            mu <- rep(mu, n)
            return(-dmvnorm(y, mu, vv, log = T))
        }
        o <- optim(foo, p = start, lower = lower, upper = upper, 
            method = "L")
            	ml.vcv <- exp(o$par[1]) * ouMatrix(vcvMat, exp(o$par[2]))
  			root.state <- phylogMean(ml.vcv, y)
  			results <- list(lnl = -o$value, root.state = root.state, beta = exp(o$par[1]), alpha = exp(o$par[2]));
  		}	
  		else
  		{
  			node=Ntip(phy) + 1
			relations_num <- monoClade(phy, node)
			originTimeFocal <- allTimes[which(names(allTimes) == node)]
			times <- c()
			for(u in 1:length(allTimes)) times[which(phy$edge[,1] == names(allTimes)[u])] <- allTimes[u]
			branches <- match(relations_num, phy$edge[,2])
			if(any(is.na(branches))) branches <- branches[complete.cases(branches)]
			start <- log(0.5)
			upper <- log(exp(1))
			lower <- log(1e-8)
			
 			foo <- function(x) {
  				tr <- OUtree(phy, alpha=exp(x), branches=branches, times=times, originTime=originTimeFocal)
				return(-phyIC(tr, y)$logLikelihood)
				}
			o <- optim(foo, p = start, lower = lower, upper = upper, method = "L");
			mlPhy <- OUtree(phy, alpha=exp(o$par), branches=branches, times=times, originTime=originTimeFocal)
			xx <- phyIC(mlPhy, y)
			results <- list(lnl = xx$logLikelihood, root.state = as.numeric(xx$phyloMean), beta = xx$sigmaSq, alpha = exp(o$par));
			if(returnPhy == TRUE) results$phy <- xx$rescaledPhylo

 		
 		}
  			
  			results$aic <- 2 * k - 2 * results$lnl
			results$aicc <- aicc(results$lnl, k, Ntip(phy))
  			results$k <- k
  			return(results)  
		}
}