#' ouMatrix function from Graham Slater
#'
#' ouMatrix function from Graham Slater
#' @param vcvMatrix
#' @param alpha
#' @keywords 
#' @return list phylogeny and matrix of parameters
#' @export
#' @examples
#' ouMatrix()

ouMatrix <- function (vcvMatrix, alpha) 
{
    vcvDiag <- diag(vcvMatrix)
    diagi <- matrix(vcvDiag, nrow = length(vcvDiag), ncol = length(vcvDiag))
    diagj <- matrix(vcvDiag, nrow = length(vcvDiag), ncol = length(vcvDiag), 
        byrow = T)
    Tij = diagi + diagj - (2 * vcvMatrix)
    vcvRescaled = (1/(2 * alpha)) * exp(-alpha * Tij) * (1 - 
        exp(-2 * alpha * vcvMatrix))
    return(vcvRescaled)
}
