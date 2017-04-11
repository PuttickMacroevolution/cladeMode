#' arbutus modified function from arbutus
#'
#' arbutus modified function from arbutus
#' @param input
#' @param testData
#' @keywords 
#' @return
#' @export
#' @examples
#' arbutus()

arbutus <- function(input, testData){
	unit.tree <- make_unit_tree.phylo(input, data=testData)
	obs <- calculate_pic_stat(unit.tree, stats=NULL)
	simdat <- simulate_char_unit(unit.tree, nsim=1000)
	sim <- calculate_pic_stat(simdat, stats=NULL)
	res <- compare_pic_stat(obs, sim)
	return(pvalue_arbutus(res))
	}
