#' Deviance residuals for NB/Poisson GLM
#'
#' @param fit edgeR glmFit object
#'
#' @return residual matrix
#' @export
devianceResiduals <- function(fit){

  y <- fit$counts
  mu <- fit$fitted.values
  theta <- 1 / fit$dispersion

  res <- matrix(0, nrow(y), ncol(y))

  for(i in seq_len(nrow(y))){

    yy <- y[i,]
    mm <- mu[i,]

    if(theta[i] == Inf){
      d <- sqrt(pmax(stats::poisson()$dev.resids(yy, mm, 1), 0))
    } else {
      d <- sqrt(pmax(MASS::negative.binomial(theta[i])$dev.resids(yy, mm, 1), 0))
    }

    res[i,] <- ifelse(yy > mm, d, -d)
  }

  res
}
