#' make_unit_tree.phylo function from arbutus
#'
#' make_unit_tree.phylo function from arbutus
#' @param x
#' @param data
#' @keywords 
#' @return
#' @export
#' @examples
#' make_unit_tree.phylo()

make_unit_tree.phylo <- function(x, data, ...) {
	require(arbutus)
	pics <- pic(data, x, var.contrasts=TRUE)
	unit.tree <- list(phy=x, data=data, pics=pics)
	class(unit.tree) <- "unit.tree"
 	unit.tree
}
