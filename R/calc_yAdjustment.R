#' Calculate y axis adjustment (two comparisons only)
#'
#' Calculates the y axis adjustment positions for placement invisible points to improve plotting space to allow the placement of brackets.
#' 
#' Obs: The groups MUST come in the order that you want them to appear in the plot.
#' @param values_df 'Melted' dataframe (long format) containing the actual values that will be plotted.
#' @param values_col Name of the column that contains the plotted values.
#' @param facet_var_value Name of the facet variable (column) in the values dataframe.
#' @param k Constant value to be added for a standard distance between the brackets. Adjust according to expression values and number of statistical comparisons.
#' @return The same values_df with an extra column called 'yadjustment' with y axis values.
#' @export
calc_yAdjustment <- function(values_df,
                             values_col,
                             facet_var_value,
                             k) {
  

  values_df$yadj <- values_df[[values_col]] + k + minMaxNorm(values_df[[values_col]], min(values_df[[values_col]]), max(values_df[[values_col]]))
  
  return(values_df)
  
}

#' Calculate y axis adjustment (>2 comparisons)
#'
#' Calculates the y axis adjustment positions for placement invisible points to improve plotting space to allow the placement of brackets.
#' 
#' Obs: The groups MUST come in the order that you want them to appear in the plot.
#' @param annotation_df Annotation dataframe to be used for plotting stats and brackets.
#' @param values_df 'Melted' dataframe (long format) containing the actual values that will be plotted.
#' @param values_col Name of the column that contains the plotted values. Default value: "Expression".
#' @param facet_var_value Name of the facet variable (column) in the values dataframe. Default value: "Gene".
#' @param facet_var_annot Name of the facet variable (column) in the annotation dataframe. Default value: "Gene".
#' @param group_var_annot Name of the grouping variable (column) in the annotation dataframe. Default value: "Group2".
#' @param base_var_annot Name of the base grouping variable (column) in the annotation dataframe. This is the group of the reference level. Default value: "Group1".
#' @param stats_col Name of the statistics variable (column) in the annotation dataframe. Default value: "padj".
#' @param k Proportion value (of the max value in plot) to be added for distance between the brackets. Adjust according to expression values and number of statistical comparisons.
#' @return The same values_df with an extra column called 'yadjustment' with y axis values.
#' @export
calc_yAdjustment2 <- function(annotation_df,
                              values_df,
                              values_col = "Expression",
                              facet_var_value = "Gene",
                              facet_var_annot = "Gene",
                              group_var_annot = "Group2",
                              base_var_annot = "Group1",
                              stats_col = "padj",
                              k) {
  
  annotation_df$yposition <- 0
  
  index <- 1:choose(length(unique(c(annotation_df[[group_var_annot]],annotation_df[[base_var_annot]]))),2)
  index <- index[1:(length(index))]
  names(index) <- unique(paste(annotation_df[[group_var_annot]], annotation_df[[base_var_annot]], sep = "_vs_"))
  
  index_list <- list()
  
  for (i in 1:length(unique(annotation_df[[facet_var_annot]]))) {
    
    facet_var <- unique(annotation_df[[facet_var_annot]])[i]
    index_list[[facet_var]] <- index
    groups <- names(index)
    
    for (comparison in groups) {
      
      group <- strsplit(comparison, "_vs_")[[1]][1]
      base <- strsplit(comparison, "_vs_")[[1]][2]
      
      stats_value <- annotation_df[annotation_df[[facet_var_annot]]==facet_var&annotation_df[[group_var_annot]]==group&annotation_df[[base_var_annot]]==base,stats_col]
      
      if ((stats_value > 0.05 | is.na(stats_value))) {
        
        index_list[[facet_var]] <- index_list[[facet_var]][names(index_list[[facet_var]]) != comparison]
        
        if (length(index_list[[facet_var]]) != 0) {
          
          index_list[[facet_var]] <- setNames(1:length(index_list[[facet_var]]), names(index_list[[facet_var]]))
          
          
        } else {
          
          index_list <- index_list[names(index_list) != facet_var]
          
        }
        
        
      }
      
    }
    
  }
  
  values_df$yadj <- 0
  
  for (facet in unique(values_df[[facet_var_annot]])) {
    max_value <- max(values_df[values_df[facet_var_value]==facet,values_col])
    values_df[values_df[[facet_var_annot]]==facet, "yadj"] <- max_value + (k * max_value * length(index_list[[facet]]))
  }
  
  return(values_df)
  
}