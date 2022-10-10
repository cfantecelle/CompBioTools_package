#' Check DESeq2 results
#'
#' Compares results from differential expression analysis and a data.frame with previously selected genes for analysis.
#' @param results_list A named list containing results objects from DE pipelines.
#' @param genes_df A dataframe of genes of interest where the column with gene symbols is named `hgnc_symbol`.
#' @param from If the results came from `'deseq2'` or `'limma'`. Default = `'deseq2'`
#' @return The function check every gene of `genes_df` against each of the results and aggregates log2 fold changes and adjusted p-values into a single data.frame for each result of the list. Results are rounded to two digits.
#' @export
geneCheck <- function(results_list,
                      genes_df,
                      from = 'deseq2') {
  
  if (from == 'deseq2') {
    
    geneCheck_results <- genes_df
    for(i in 1:length(results_list)) {
      if(class(results_list[[i]]) == "DESeqResults") {
        col_fc <- paste0(names(results_list[i]), "_log2FC")
        col_padj <- paste0(names(results_list[i]), "_padj")
        for(g in 1:length(genes_df$hgnc_symbol)) {
          if(genes_df$hgnc_symbol[g] %in% rownames(results_list[[i]])) {
            geneCheck_results[g,col_fc] <- round(results_list[[i]][genes_df$hgnc_symbol[g],]$log2FoldChange, 2)
            geneCheck_results[g,col_padj] <- as.numeric(formatC(results_list[[i]][genes_df$hgnc_symbol[g],]$padj, format = "e", digits = 2))
          } else {
            geneCheck_results[g,col_fc] <- NA
            geneCheck_results[g,col_padj] <- NA
          }
        }
      } else {
        message(paste0("The result ", names(geneCheck_results[i], "is not a DESeqResults object.")))
      }
    }
    return(geneCheck_results)
    
  } else if (from == 'limma') {
    
    geneCheck_results <- genes_df
    
    geneCheck_results <- genes_df
    for(i in 1:length(results_list)) {
      if(class(results_list[[i]]) == "data.frame") {
        col_fc <- paste0(names(results_list[i]), "_log2FC")
        col_padj <- paste0(names(results_list[i]), "_padj")
        for(g in 1:length(genes_df$hgnc_symbol)) {
          if(genes_df$hgnc_symbol[g] %in% rownames(results_list[[i]])) {
            geneCheck_results[g,col_fc] <- round(results_list[[i]][genes_df$hgnc_symbol[g],]$logFC, 2)
            geneCheck_results[g,col_padj] <- as.numeric(formatC(results_list[[i]][genes_df$hgnc_symbol[g],]$adj.P.Val, format = "e", digits = 2))
          } else {
            geneCheck_results[g,col_fc] <- NA
            geneCheck_results[g,col_padj] <- NA
          }
        }
      } else {
        message(paste0("The result ", names(geneCheck_results[i], "is not a Limma results object (data.frame class).")))
      }
    }
    return(geneCheck_results)
    
  } else {
    stop('Results format not yet accepted.')
  }
  
}