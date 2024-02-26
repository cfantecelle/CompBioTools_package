#' Calculate y axis position (two comparisons only)
#'
#' Calculates the y axis position for placement of brackets when adding stats to faceted plots from a melted matrix.
#' 
#' Obs: The groups MUST come in the order that you want them to appear in the plot.
#' @param annotation_df Annotation dataframe to be used for plotting stats and brackets.
#' @param values_df 'Melted' dataframe (long format) containing the actual values that will be plotted.
#' @param values_col Name of the column that contains the plotted values.
#' @param facet_var_value Name of the facet variable (column) in the values dataframe.
#' @param facet_var_annot Name of the facet variable (column) in the annotation dataframe.
#' @param group_var_annot Name of the grouping variable (column) in the annotation dataframe.
#' @param k Proportion value (of the max value in plot) to be added for distance between the brackets. Adjust according to expression values and number of statistical comparisons.
#' @return The same annotation_df with an extra column called 'yposition' with y axis values.
#' @export
calc_yPosition <- function(annotation_df,
                           values_df,
                           values_col,
                           facet_var_value,
                           facet_var_annot,
                           group_var_annot, 
                           k) {
  
  annotation_df$yposition <- 0
  
  index <- 1:length(unique(annotation_df[[group_var_annot]]))
  names(index) <- unique(annotation_df[[group_var_annot]])
  message("The order of the groupings is as follows:\n")
  print(index)
  
  for (i in 1:length(annotation_df[[facet_var_annot]])) {
    group <- annotation_df[[group_var_annot]][i]
    max_value <- max(values_df[[values_col]][values_df[[facet_var_value]]==annotation_df[[facet_var_annot]][i]])
    annotation_df$yposition[i] <- max_value + (k * max_value * index[[group]])
  }
  
  return(annotation_df)
  
}



#' Calculate y axis position (>2 comparisons)
#'
#' Calculates the y axis position for placement of brackets when adding stats to faceted plots from a melted matrix.
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
#' @param alpha Significance threshold.
#' @param k Proportion value (of the max value in plot) to be added for distance between the brackets. Adjust according to expression values and number of statistical comparisons.
#' @param decrease Integer number to subtract from total number of comparisons. This is for use in case you are plotting stats for less than the maximum amount of comparisons between the groups. Default is `0`.
#' @return The same annotation_df with an extra column called 'yposition' with y axis values.
#' @export
calc_yPosition2 <- function(annotation_df,
                            values_df,
                            values_col = "Expression",
                            facet_var_value = "Gene",
                            facet_var_annot = "Gene",
                            group_var_annot = "Group2",
                            base_var_annot = "Group1",
                            stats_col = "padj",
                            alpha = 0.05,
                            k,
                            decrease = 0) {
  
  annotation_df$yposition <- 0
  
  if (length(unique(annotation_df[[base_var_annot]])) > 1) {
    
    index <- 1:choose(length(unique(c(annotation_df[[group_var_annot]],annotation_df[[base_var_annot]]))),2)
    index <- index[1:(length(index)-decrease)]
    names(index) <- unique(paste(annotation_df[[group_var_annot]], annotation_df[[base_var_annot]], sep = "_vs_"))
    message("The order of the groupings is as follows:\n")
    print(index)
    
  } else {
    
    index <- 1:length(unique(annotation_df[[group_var_annot]]))
    names(index) <- unique(paste(annotation_df[[group_var_annot]], annotation_df[[base_var_annot]], sep = "_vs_"))
    message("The order of the groupings is as follows:\n")
    print(index)
    
  }
  
  index_list <- list()
  
  for (i in 1:length(unique(annotation_df[[facet_var_annot]]))) {
    
    facet_var <- unique(annotation_df[[facet_var_annot]])[i]
    index_list[[facet_var]] <- index
    groups <- names(index)
    
    for (comparison in groups) {
      
      group <- strsplit(comparison, "_vs_")[[1]][1]
      base <- strsplit(comparison, "_vs_")[[1]][2]
      
      stats_value <- annotation_df[annotation_df[[facet_var_annot]]==facet_var&annotation_df[[group_var_annot]]==group&annotation_df[[base_var_annot]]==base,stats_col]
      
      if (stats_value > 0.05) {
        
        index_list[[facet_var]] <- index_list[[facet_var]][names(index_list[[facet_var]]) != comparison]
        
        if (length(index_list[[facet_var]]) != 0) {
          
          index_list[[facet_var]] <- setNames(1:length(index_list[[facet_var]]), names(index_list[[facet_var]]))
          
          
        } else {
          
          index_list <- index_list[names(index_list) != facet_var]
          
        }
        
        
      }
      
    }
    
  }
  
  
  for (i in 1:length(annotation_df[[facet_var_annot]])) {
    
    # Extracting group comparison and facet ver
    group1 <- annotation_df[i, base_var_annot]
    group2 <- annotation_df[i, group_var_annot]
    group <- paste(group2, group1, sep = "_vs_")
    
    facet_var <- annotation_df[i, facet_var_annot]
    
    # Calculating y position of brackets
    if (length(index_list[[facet_var]][group]) > 0) {
      
      max_value <- max(values_df[values_df[facet_var_value]==facet_var,values_col])
      
      
      
      if (annotation_df[i, stats_col] < alpha) {
        
        annotation_df$yposition[i] <- max_value + (k * max_value * index_list[[facet_var]][group])
      
      } else {
        
        annotation_df$yposition[i] <- NA
        
        
        }
      
      
    } else {
      
      annotation_df$yposition[i] <- NA
      
    }
    
  }
  
  values_df$yadj <- 0
  
  for (facet in unique(values_df[[facet_var_annot]])) {
    max_value <- max(values_df[values_df[facet_var_value]==facet_var,values_col])
    values_df[values_df[[facet_var_annot]]==facet, "yadj"] <- max_value + (k * max_value * length(index_list[[facet]]))
  }
  
  return(annotation_df)
  
}