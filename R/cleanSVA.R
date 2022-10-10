#' Regress out batch/unwanted variation using SVA
#'
#' @description 
#' Remove batch/unwanted effects from count data.
#' 
#' @details
#' Takes the result from Surrogate Variable Analysis (see [`svaseq()`]) and use it to regress out
#' the unwanted effects from count data. 
#' This method was used by Jaffe and coleagues in the following paper: [10.1186/s12859-015-0808-5](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0808-5).
#' 
#' @param y The counts matrix from which to remove the effects.
#' @param mod Model matrix to be tested (full model of the known effects). Must be the same used at [`svaseq()`]
#' @param svaobj The result from [`svaseq()`].
#' @returns 'Clean' matrix with unwanted effects removed.
#' @export
cleanSVA <- function(y, mod, svaobj) {
  
  X <- cbind(mod, svaobj$sv)
  Hat <- solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(y))
  P <- ncol(mod)
  
  cleanMatrix <- y-t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P),])
  
  return(cleanMatrix)
}