#' Run full CTRLseq pipeline
#'
#' @param counts gene x sample count matrix
#' @param group factor (control/case)
#' @param n_factor number of latent factors
#' @param beta KL divergence weight
#' @param seed random seed
#'
#' @return DEG result table
#' @export

if (requireNamespace("tensorflow", quietly = TRUE)) {
  `%as%` <- tensorflow::`%as%`
}

run_CTRLseq <- function(counts, group, n_factor = 2, beta = 0.1, seed = 123){

  set.seed(seed)
  Sys.setenv(TF_ENABLE_ONEDNN_OPTS = 0)

  if (!requireNamespace("tensorflow", quietly = TRUE)) {
    stop("tensorflow required")
  }

  tensorflow::tf$random$set_seed(seed)
  reticulate::py_run_string(paste0("import random; random.seed(", seed, ")"))

  library(edgeR)
  library(keras3)
  library(tensorflow)

  `%as%` <- tensorflow::`%as%`

  obj <- CTRLseq(
    counts = counts,
    group = group,
    n.factor = n_factor,
    beta = beta,
    seed = seed
  )

  res <- DEGs_CTRL(obj)

  return(res)
}
