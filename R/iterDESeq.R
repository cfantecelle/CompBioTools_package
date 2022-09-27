#' DESeq with maxit argument
#'
#' @description 
#' Runs DESeq2 with the different maximum iterations value.
#' 
#' @details
#' This function executes the same steps executed by the `DESeq2()` wrapper function, allowing you to change the maximum iterations value. These steps are:
#' 
#' `dds <- estimateSizeFactors(dds)`
#' `dds <- estimateDispersions(dds)`
#' `dds <- nbinomWaldTest(dds, maxit = n)`
#' 
#' Where `n` will be the chosen maximum iterations number in this function's `maxit` parameter.
#' 
#' @param obj The DESeq object on which to execute DESeq2()
#' @param maxit Number of maximum iterations to ther nbinomWaldTest test in DESeq2. Default is 500.
#' @returns The same DESeq object after running all the functions above.
#' @export
iterDESeq <- function(obj, maxit = 500) {
  
  obj <- DESeq2::estimateSizeFactors(obj)
  obj <- DESeq2::estimateDispersions(obj)
  obj <- DESeq2::nbinomWaldTest(obj, maxit = maxit)
  
  return(obj)
}