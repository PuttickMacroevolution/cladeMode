#' split.vcv.node
#'
#' split.vcv.node produces two VCV matrices that represent a monophyletic subclade and the remaining phylogeny respectively
#' @param phy input phylogeny
#' @param node 
#' @keywords 
#' @return vcvMatrix vcv matrix
#' @export
#' @examples
#' set.seed(30)
#' #simulate tree with 20 species
#' simTree <- rcoal(20)
#' split.vcv.node(simTree, 22)

split.vcv.node <- function(phy, node){
  mat <- vcv(phy)
  forty <- extract.clade(phy, node)
  vcv2 <- vcv(forty)
  relations <- monoClade(phy=phy, node)
  rel <- relations[which(relations <= Ntip(phy))]
  nodesRows <- match(phy$tip.label[rel], rownames(mat))
  baseMat <- mat
  reArr <-match(rownames(mat), rownames(vcv2))
  ff <- vcv(forty)[reArr, reArr]
  ff[which(is.na(ff) == T)] <- 0
  matOne <- baseMat - ff
  vcv2 <- mat
  vcv2[,] <- 0
  matTwo <- vcv2 + ff
  mats <- list(matOne, matTwo)
  names(mats) <- c("base", node)
  vcvMatrix <- mats
  return(vcvMatrix)
}
