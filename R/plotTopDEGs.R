#' Plot top differentially expressed genes
#'
#' Takes a list o results from differential expression pipelines and plot heatmaps of the top `n` differentially expressed
#' genes (DEGs) based on fold change value and plots them using [ComplexHeatmap] into a single heatmap or list of heatmaps.
#' The gene list can be calculated using positive fold changes, negative fold changes or absolute values.
#' @param results A data.frame or list of data.frames containing differential expression results.
#' @param counts_matrix Normalized and log transformed matrix of counts for plotting. Should be already centered.
#' @param n_genes Number of top DEGs to select. Default = `20`.
#' @param which A string indicating which direction to choose, `'up'` for overrepresented genes, `'down'` for downrepresented genes and `'abs'` for absolute values (both ways). Default: `'abs'`.
#' @param logFC_col The name of the column where the fold changes values are stored. Not a string.
#' @param pval_col Name of the column where the p-values (adjusted) are stored. Not a string.
#' @param groups A `data.frame` with results names as rownames (the same as names of results in list), and two columns in this order: one indicating the baseline group and one for the compared group. These must be the sames as in `pdata` argument `data.frame`.
#' @param pdata A `data.frame` containing `sample` and `groups` columns.
#' @param sample_col Name of column with sample information. Not a string.
#' @param group_col Name of column with group information. Not a string.
#' @param split_col `TRUE` or `FALSE`. Specifies wether to split the columns of the heatmaps based on sample information.
#' @return A Heatmap or List of Heatmaps as per the [ComplexHeatmap] package.
#' @export
plotTopDEGs <- function(results,
                        counts_matrix,
                        n_genes = 20,
                        which = 'abs',
                        logFC_col,
                        pval_col,
                        groups,
                        pdata,
                        sample_col = sample,
                        group_col = group,
                        split_col = TRUE) {
  

  
  ## Checks
  
  # Results df
  if (class(results) != 'data.frame' & class(results) != 'DESeqResults') {
    if (class(results) == 'list') {
      if (length(unique(unlist(lapply(results, class)))) > 1) {
        return(stop("Results should all be the same type."))
        print(results)
      } else if (unique(unlist(lapply(results, class))) != 'data.frame' &
                 unique(unlist(lapply(results, class))) != 'DESeqResults') {
        return(stop("List provided is of neither data.frame or DESeqResults types."))
      }
    } else {
      return(stop("Results provided are not either a data.frame, a DESeqResults object or a list of those."))
      print(results)
    }
  }
  
  # Reordering pdata de matrix
  counts_matrix <- counts_matrix[,order(colnames(counts_matrix))]
  pdata <- pdata %>% arrange({{sample_col}})
  rownames(pdata) <- pdata %>% pull({{sample_col}})
  
  if (!identical(rownames(pdata), colnames(counts_matrix))) {
    return(stop("Column names of counts_matrix do not match the samples names from pdata."))
  }
  
  if (class(results) == 'list') {
    
    # Ensuring only significant values
    res_signif <- results %>% map(filter, {{pval_col}} < 0.05)
    
    # Sorting by log2FC
    if (which == 'up') {
      res_topN <- res_signif %>% map(arrange, desc({{logFC_col}}))
    } else if (which == 'down') {
      res_topN <- res_signif %>% map(arrange, {{logFC_col}})
    } else if (which == 'abs') {
      res_topN <- res_signif %>% map(arrange, desc(abs({{logFC_col}})))
    } else {
      return(stop(glue("'{which}' is an invalid option for argument 'which'.")))
    }
    
    # Chopping
    res_topN <- res_topN %>% map(function(x) { x[1:n_genes, , drop = F]})
    
    # Subsetting matrix
    matrices <- list()
    
    for (res in names(results)) {
      comps <- as.vector(t(groups[res,]))
      samples <- pdata %>% filter({{group_col}} %in% comps) %>% pull({{sample_col}})
      temp <- counts_matrix[rownames(res_topN[[res]]),samples]
      matrices[[res]] <- temp
      
    }
    
    # Creating heatmaps
    hm_list <- HeatmapList()
    
    for (res in names(results)) {
      
      # Creating col_function
      col_fun <- colorRamp2(c(min(counts_matrix), 0, abs(min(counts_matrix))),
                            c("#00c7ff", "#000000", "#fff300"))
      
      # Creating annotation df
      annot_df <- pdata %>% filter({{sample_col}} %in% colnames(matrices[[res]])) %>% select({{group_col}})
      colnames(annot_df) <- "Group"
      
      # Creating annotation
      top_annot_horiz <- HeatmapAnnotation(df = annot_df, show_annotation_name = FALSE,
                                           simple_anno_size = unit(.4, "cm"),
                                           annotation_legend_param = list(nrow = 2, title_position = "topcenter"),
                                           annotation_name_align = TRUE,
                                           show_legend = TRUE)
      
      hm_temp <- Heatmap(matrices[[res]], name = res,
                         column_title = "",
                         column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                         row_names_gp = gpar(fontsize = 10, fontface = "italic"),
                         col = col_fun,
                         top_annotation = top_annot_horiz,
                         column_split = annot_df$Group,
                         show_row_names = TRUE,
                         show_column_names = FALSE,
                         show_row_dend = FALSE,
                         show_column_dend = FALSE,
                         cluster_column_slices = FALSE,
                         show_parent_dend_line = FALSE,
                         use_raster = FALSE,
                         cluster_rows = TRUE,
                         cluster_columns = TRUE,
                         clustering_distance_columns = "euclidean",
                         heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"))
      
      hm_list <- hm_list + hm_temp
      
    }
    
    return(hm_list)
    
  } else {
    
    message("No method for single objects yet...")
    
  }
  
}
  
