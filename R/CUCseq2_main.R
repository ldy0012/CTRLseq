#' CUCseq2 main model
#'
#' @param counts gene x sample count matrix
#' @param group factor (control/case)
#' @param n.factor number of latent factors
#' @param beta KL divergence weight
#' @param seed random seed
#'
#' @return list containing fitted model, latent factors, and DGE object
#' @export
CUCseq2 <- function(counts, group, n.factor = 2, beta = 0.1, seed = 123){

  design <- model.matrix(~group)

  dge <- edgeR::DGEList(counts = counts)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)

  fit <- edgeR::glmFit(dge, design)

  E <- residuals(fit, type = "deviance")

  F <- estimate_latent_factors(E, n.factor, beta = beta, seed = seed)
  F <- scale(F)

  design2 <- cbind(design, F)

  dge2 <- edgeR::DGEList(counts = counts)
  dge2 <- edgeR::calcNormFactors(dge2)
  dge2 <- edgeR::estimateDisp(dge2, design2)

  fit2 <- edgeR::glmFit(dge2, design2)

  list(
    fit = fit2,
    factors = F,
    dge = dge2
  )
}
