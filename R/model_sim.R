#' model_sim
#'
#' summarise simulations
#' @param phy input phylogeny
#' @param sim_data
#' @param n_core default=1
#' @keywords 
#' @return
#' @export
#' @examples
#' model_sim()

model_sim <- function(phy, sim_data, n_core=1, nReps=100) {
	bm_res <- mclapply(1:nReps, mc.cores=n_core, function(x) cladeMode(phy, sim_data[,x], mode="BM", contrasts=T))
	bm_res <- matrix(unlist(bm_res), ncol=nReps)
	eb_res <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], mode="EB", contrasts=T))
	eb_res <- matrix(unlist(eb_res), ncol=nReps)
	ou_res <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], mode="OU", contrasts=T))
	ou_res <- matrix(unlist(ou_res), ncol=nReps)
	minSize_0.05 <- ceiling(Ntip(phy) * 0.05)
	nested_eb_res <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], minSize=minSize_0.05, mode="nestedEB", contrasts=T))
	nested_eb_res <- matrix(unlist(nested_eb_res), ncol=nReps)
	nested_eb_rate_res <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], minSize=minSize_0.05, mode="nestedEBRate", contrasts=T))
	nested_eb_rate_res <- matrix(unlist(nested_eb_rate_res), ncol=nReps)
	minSize_0.25 <- ceiling(Ntip(phy) * 0.25)
	nested_ou_res_0.25 <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], minSize=minSize_0.25, mode="nestedOU", contrasts=T))		
	nested_ou_res_0.25 <- matrix(unlist(nested_ou_res_0.25), ncol=nReps)
	nested_eb_res_0.25 <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], minSize=minSize_0.25,mode="nestedEB", contrasts=T))
	nested_eb_res_0.25 <- matrix(unlist(nested_eb_res_0.25), ncol=nReps)
	nested_eb_rate_res_0.25 <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], minSize=minSize_0.25,mode="nestedEBRate", contrasts=T))
	nested_eb_rate_res_0.25 <- matrix(unlist(nested_eb_rate_res_0.25), ncol=nReps)
	nested_rate_res_0.25 <- mclapply(1:nReps, mc.cores=n_core,function(x) cladeMode(phy, sim_data[,x], minSize=minSize_0.25, mode="rateShift", contrasts=T))
	nested_rate_res_0.25 <- matrix(unlist(nested_rate_res_0.25), ncol=nReps)
	final_res <- list(bm_res, eb_res, ou_res, nested_eb_res, nested_eb_rate_res, 
	nested_ou_res_0.25, nested_eb_res_0.25, nested_eb_rate_res_0.25, nested_rate_res_0.25)
	return(final_res)
}