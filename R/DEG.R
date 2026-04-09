#' Extract DEGs from CUCseq2 result
#'
#' @param obj output from CUCseq2
#'
#' @return DEG result table
#' @export
DEGs_CUC <- function(obj){

  lrt <- edgeR::glmLRT(obj$fit, coef = 2)
  tab <- edgeR::topTags(lrt, n = nrow(obj$fit$counts))$table

  tab
}
