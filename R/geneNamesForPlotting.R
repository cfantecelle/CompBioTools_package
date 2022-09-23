#' Gene names for plotting
#'
#' Create labels for plotting to improve understanding of data by allowing the readers to identify genes by their common aliases.
#' @param genes_df A dataframe containing the genes analysed. Requires 'name' column with common name and 'hgnc_symbol' with the official gene symbol.
#' @return The same dataframe with an extra column (names_plot column) with names in the format: SYMBOL(Alias).
#' @export
geneNamesForPlotting <- function(genes_df) {
  
  # Takes dataframe with genes of interest and make a column with names for plotting
  # Requires 'name' column with common name and 'hgnc_symbol' with the official gene symbol
  
  genes_df$names_plot <- ""
  
  for (i in 1:length(genes_df$hgnc_symbol)) {
    if (genes_df$hgnc_symbol[i] != genes_df$name[i]) {
      genes_df$names_plot[i] <- paste(genes_df$hgnc_symbol[i], " (", genes_df$name[i], ")", sep = "")
    } else {
      genes_df$names_plot[i] <- genes_df$hgnc_symbol[i]
    }
  }
  
  return(genes_df)
}