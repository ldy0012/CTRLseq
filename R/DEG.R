#' Extract DEGs from CTRLseq result
#'
#' @param obj output from CTRLseq
#'
#' @return DEG result table
#' @export
DEGs_CTRL <- function(obj){

  lrt <- edgeR::glmLRT(obj$fit, coef = 2)
  tab <- edgeR::topTags(lrt, n = nrow(obj$fit$counts))$table

  tab
}
