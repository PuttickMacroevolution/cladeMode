#' independent contrasts calculation
#'
#' independent contrasts calculation
#' @param phy
#' @param yDat
#' @keywords 
#' @return list of parameters
#' @export
#' @examples
#' phyIC()

phyIC <- function(phy, yDat) {
	

		bL <- phy$edge.length
		tips <- which(phy$edge[,2] <= Ntip(phy))
		tipsMat <- phy$edge[tips,]
		y <- trait <- yDat
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
	
		
		recons <- rep(NA, length(c(Ntip(phy)+1, phy$edge[, 2])))
		names(recons) <- c(Ntip(phy)+1, phy$edge[, 2])
		correctOrder <- as.numeric(names(recons)[which(as.numeric(names(recons)) <= Ntip(phy))])
		reArr <- order(correctOrder)
		recons[which(as.numeric(names(recons)) <= Ntip(phy))] <- y[reArr]

		count <- 1
		
		while(1 == 1) {
		stillToDo <- which(is.na(recons))
		if(length(stillToDo) == 0) break()
		incom <- as.numeric(names(recons)[stillToDo])
		incom <- rev(incom)
		allIn <- which(phy$edge[,1] == incom[count])
		doneTest <- match(phy$edge[allIn, 2], names(recons))
		decCom <- any(is.na(recons[doneTest]))	
		if(decCom == F) {
			yx <- recons[doneTest]
			correct <- rescaled[allIn]
			cor <- prod(correct)/ correct
			rec <- weighted.mean(yx, cor)
			insert <- match(incom[count], names(recons))
			recons[insert] <- rec
			count <- 1
			} else {
			count <- count + 1
			}
		}
	
		ints <- (Ntip(phy) + 1) : (Ntip(phy) + Nnode(phy))
		yx <- recons[-1]
		picCont <- sapply(ints, function(p) {
			partners <- which(phy$edge[,1] == p)
			contrasts <- diff(rev(yx[partners])) 
			scaledContrasts <- contrasts / (sqrt(sum(rescaled[partners])))
			variance <- sum(rescaled[partners])
			c(contrasts, scaledContrasts, variance)
		})
		
		part.root <- partners <- which(phy$edge[,1] == ints[1])
		rootVar <- prod(rescaled[partners])/(sum(rescaled[partners]))
	
		picCont <- as.matrix(t(picCont))
		colnames(picCont) <- c("contrasts", "scaledContrasts", "variance")
		rownames(picCont) <- ints
	
		vars <- picCont[,3]
		rawCon <- picCont[,1]
		vars <- c(rootVar, vars)
		sigmaSq <- sum(c(0, rawCon)^2 / vars) / n 
		phyloMean <- recons[1]
		
		likLog <- -0.5 * (n * log(2*pi*sigmaSq) + sum(log(vars) + (c(0, picCont[,1])^2 / (sigmaSq * vars))))	
	
		rescaledPhylo <- phy
		rescaledPhylo$edge.length <- rescaled
	
		return(list(contrasts = picCont, sigmaSq = sigmaSq, phyloMean = phyloMean, logLikelihood = likLog, rescaledPhylo=rescaledPhylo, ancestralStates = recons, v0 = rootVar))
	}
