#' Run full CUCseq2 pipeline
#'
#' @param counts gene x sample count matrix
#' @param group factor (control/case)
#' @param n_factor number of latent factors
#' @param beta KL divergence weight
#' @param seed random seed
#'
#' @return DEG result table
#' @export
run_CUCseq2 <- function(counts, group, n_factor = 2, beta = 0.1, seed = 123){

  obj <- CUCseq2(
    counts = counts,
    group = group,
    n.factor = n_factor,
    beta = beta,
    seed = seed
  )

  res <- DEGs_CUC(obj)

  return(res)
}
