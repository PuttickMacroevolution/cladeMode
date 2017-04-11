#' ACDC.prior function from Geiger
#'
#' ACDC.prior function from Geiger
#' @param phy
#' @keywords 
#' @return
#' @export
#' @examples
#' ACDC.prior()

ACDC.prior <- function (phy, decrease.max = 1e-02, increase.max = 1e+02) 
{
    if (is.ultrametric(phy)) {
        max.bt <- max(branching.times(phy))
    }
    else {
        max.bt <- max(BranchingTimesFossil(phy))
    }
    prior.min <- log(decrease.max)/max.bt
    prior.max <- log(increase.max)/max.bt
    return(list(min = prior.min, max = prior.max))
}
