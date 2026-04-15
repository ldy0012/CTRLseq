#' CTRLseq main model
#'
#' @param counts gene x sample count matrix
#' @param group factor (control/case)
#' @param n.factor number of latent factors
#' @param beta KL divergence weight
#' @param seed random seed
#'
#' @return list containing fitted model, latent factors, and DGE object
#' @export
CTRLseq <- function(counts, group, n.factor = 2, beta = 0.1, seed = 123){

  design <- model.matrix(~group)

  dge <- edgeR::DGEList(counts = counts, group = group)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)

  fit <- edgeR::glmFit(dge, design)

  E <- residuals(fit, type = "deviance")

  if (is.null(E)) {
    stop("Residuals are NULL. Model fitting failed.")
  }

  if (is.null(dim(E))) {
    E <- matrix(E, nrow = nrow(counts))
  }

  if (!is.matrix(E)) {
    stop("Residual matrix E is not valid.")
  }

  F <- estimate_latent_factors(E, n.factor, beta = beta, seed = seed)
  F <- scale(F)

  design2 <- cbind(design, F)

  dge2 <- edgeR::DGEList(counts = counts, group = group)
  dge2 <- edgeR::calcNormFactors(dge2)
  dge2 <- edgeR::estimateDisp(dge2, design2)

  fit2 <- edgeR::glmFit(dge2, design2)

  list(
    fit = fit2,
    factors = F,
    dge = dge2
  )
}
