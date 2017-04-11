#' phylogMean function from Graham Slater
#'
#' phylogMean function from Graham Slater
#' @param phyvcv
#' @param data
#' @keywords 
#' @return list phylogeny and matrix of parameters
#' @export
#' @examples
#' phylogMean()

phylogMean <- function (phyvcv, data) 
{
    o <- rep(1, length(data))
    ci <- solve(phyvcv)
    m1 <- solve(t(o) %*% ci %*% o)
    m2 <- t(o) %*% ci %*% data
    return(m1 %*% m2)
}
