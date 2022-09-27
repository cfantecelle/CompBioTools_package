#' DESeq2 with maxiter argument
#'
#' Runs DESeq2 with the different maximum iterations value. To do that, it runs separately the following functions:
#' 
#' dds <- estimateSizeFactors(dds)
#' dds <- estimateDispersions(dds)
#' dds <- nbinomWaldTest(dds, maxit = n)
#' 
#' Where n will be the chosen maximum iterations number.
#' 
#' @param obj The DESeq object on which to execute DESeq2()
#' @param maxit Number of maximum iterations to ther nbinomWaldTest test in DESeq2. Default is 500.
#' @return DESeq2 object after running all the functions.
#' @export
iterDESeq2 <- function(obj, maxit = 500) {
  
  obj <- DESeq2::estimateSizeFactors(obj)
  obj <- DESeq2::estimateDispersions(obj)
  obj <- DESeq2::nbinomWaldTest(obj, maxit = maxit)
  
  return(obj)
}