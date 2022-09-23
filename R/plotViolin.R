#' Plot preliminary Violin plots of Gene Expression data.
#'
#' This function takes a list a genes and extracts values from a count matrix to produce Violin plots of the selected genes to be used as a preliminary analysis.
#' Currently supports only comparisons between two groups.
#' @param counts_matrix A count matrix with official gene symbols as rownames.
#' @param results_df A results dataframe to obtain statists results from.
#' @param genes_df A dataframe containing list of genes to be analysed. Must contain at least one column named 'hgnc_symbol' with official gene symbols. The 'name' column is optional and is used to plot official names along with their aliases in a SYMBOL(Alias) format.
#' @param results_from Where the DGE results were calculated. Currently, the possible values are 'deseq2' (default) or 'limma'.
#' @param pheno_data Dataframe containin phenotypic data about the samples analysed.
#' @param control_group Name of the baseline group for comparison in pheno_data.
#' @param control_color Color for plotting the baseline group.
#' @param test_group Name of the tested group for comparison in pheno_data.
#' @param test_color Color for plotting the tested group.
#' @param group_var Name of the group variable in pheno_data.
#' @param plots_order If specified, defines the order in which the genes must be plotted.
#' @param n_rows If specified, defines the number of rows of plots in the final plotting area.
#' @param n_cols If specified, defines the number of columns of plots in the final plotting area.
#' @param alpha Significance threshold (Default = 0.05).
#' @return A list object with three objects:
#' @return  • data = The genes_df dataframe with new gene information.
#' @return  • plot_annotation = The annotation_df used for plotting.
#' @return  • violin_plot = The violin plots generated from the data as faceted ggplot2 plots.
#' @export
plotViolin <- function(counts_matrix,
                       results_df,
                       genes_df,
                       results_from = "deseq2",
                       pheno_data,
                       sample_col = "sample",
                       control_group,
                       control_color = "#33C4FF",
                       test_group,
                       test_color = "#FD4646",
                       group_var,
                       plot_order = NULL,
                       n_rows = NULL,
                       n_cols = NULL,
                       alpha = 0.05) {
  
  # This function takes a list a genes and extracts values from a
  # count matrix to produce Violin plots of the selected genes.
  #
  # counts_matrix = a count matrix with gene names as row names.
  # results_df = results data.frame to obtain statistics.
  # results_from = where the results came from, either 'limma' or 'deseq2' (default value).
  # genes_df = a data.frame, must contain at least one column:
  #               
  #               hgnc_symbol -> official gene symbol
  #               name -> common gene name, alias (optional)
  #
  # control_group = base group the results use for comparison. should be the same value as in pheno_data group var.
  # test_group = tested group of interest. should be the same value as in pheno_data group var.
  # group_var = name of column used to store grouping variables mentioned above.
  # plots_order = array expliciting the order of genes for plotting.
  # n_rows = number of rows for plotting (if more than 1 gene).
  # n_cols = number of columns for plotting (same as above).
  # alpha = significance level for DEG definition.
  
  if (results_from == 'deseq2') {
    logFC_col <- 'log2FoldChange'
    meanExp_col <- 'baseMean'
    logFCSE_col <- 'lfcSE'
    pValue_col <- 'pvalue'
    pAdj_col <- 'padj'
  } else if (results_from == 'limma') {
    logFC_col <- 'logFC'
    meanExp_col <- 'AveExpr'
    logFCSE_col <- NULL
    pValue_col <- 'P.Value'
    pAdj_col <- 'adj.P.Val'
  } else {
    return(stop(paste0("Results data frame format ", results_from, " currently not supported or invalid.")))
  }
  
  ## Preparing data
  
  # Subsetting original matrix by gene
  genes_matrix <- counts_matrix[rownames(counts_matrix) %in% genes_df$hgnc_symbol,]
  
  # Possible samples grouping
  possible_groups = c(control_group, test_group)
  
  # Subsetting original pheno data
  pheno_sub <- pheno_data[pheno_data[,group_var] %in% possible_groups,]
  
  # Subsetting original matrix for samples of tested groups
  genes_matrix <- genes_matrix[,colnames(genes_matrix) %in% pheno_sub[,sample_col]]
  
  # Reordering both sets
  genes_df <- genes_df[order(genes_df$hgnc_symbol),]
  genes_matrix <- genes_matrix[order(rownames(genes_matrix)),]
  
  # Keeping only genes that appear in count matrix
  discarded_genes <- genes_df[genes_df$hgnc_symbol %notin% rownames(genes_matrix),]
  discarded_genes <- toString(discarded_genes$hgnc_symbol)
  genes_df <- genes_df[genes_df$hgnc_symbol %in% rownames(genes_matrix),]
  
  if (is.null(discarded_genes)) {
    message("All genes in genes_df present in counts_matrix. Proceeding.")
  } else {
    message(paste0("Only genes present in counts_matrix will be kept in the analysis. Genes removed from the analysis: ", discarded_genes, ".\n"))
  }
  
  # Adjusting names for plotting
  if ("name" %in% colnames(genes_df)) {
    genes_df$names_plot <- ""
    for (i in 1:length(genes_df$hgnc_symbol)) {
      if (genes_df$hgnc_symbol[i] != genes_df$name[i]) {
        genes_df$names_plot[i] <- paste(genes_df$hgnc_symbol[i], " (", genes_df$name[i], ")", sep = "")
      } else {
        genes_df$names_plot[i] <- genes_df$hgnc_symbol[i]
      }
    }
    # Replacing
    if (identical(genes_df$hgnc_symbol, rownames(genes_matrix))){
      rownames(genes_matrix) <- genes_df$names_plot
    } else {
      return(stop("Gene symbols in genes_df do not match with rownames of counts_matrix."))
    }
  } else {
    message("Gene names or aliases (\'name\' column) are not available. Results will be displayed using official gene symbols.")
    genes_df$names_plot <- genes_df$hgnc_symbol
  }
  
  ## Summarizing stats for genes for ease of use
  
  # Creating extra columns
  genes_df$log2fc <- 0
  genes_df$padj <- 0
  genes_df$signif <- ""
  genes_df$max <- 0
  genes_df$min <- 0
  
  # Filling values and creating significance codes
  for (i in 1:dim(genes_df)[1]) {
    gene <- genes_df$hgnc_symbol[i]
    genes_df$log2fc[i] <- results_df[gene, logFC_col]
    genes_df$padj[i] <- results_df[gene, pAdj_col]
    if (!is.na(genes_df$padj[i])) {
      if (genes_df$padj[i] < 0.05 & genes_df$padj[i] >= 0.01) {
        genes_df$signif[i] <- "*"
      } else if (genes_df$padj[i] < 0.01 & genes_df$padj[i] >= 0.001) {
        genes_df$signif[i] <- "**"
      } else if (genes_df$padj[i] < 0.001 & genes_df$padj[i] >= 0.0001) {
        genes_df$signif[i] <- "***"
      } else if (genes_df$padj[i] < 0.0001 & genes_df$padj[i] > 0) {
        genes_df$signif[i] <- "****"
      } else {
        genes_df$signif[i] <- "ns"
      }
    } else {
      genes_df$signif[i] <- "NA"
    }
    gene <- genes_df$names_plot[i]
    genes_df$max[i] <- max(genes_matrix[gene,])
    genes_df$min[i] <- min(genes_matrix[gene,])
    
  } 
  
  # Transposing matrix for melting
  counts_df <- as.data.frame(t(genes_matrix))
  
  # Check if samples match with phenotypic data
  if (identical(rownames(counts_df), pheno_sub[,sample_col])) {
    counts_df <- cbind(counts_df, pheno_sub)
  } else {
    stop("Sample names of counts_matrix and pheno_sub do not match!")
  }
  
  # Melting df for plotting
  countsMelted_df <- melt(counts_df, 
                          id.vars = colnames(pheno_sub), 
                          value.name = "Expression", 
                          variable.name = "Gene")
  
  # Creating annotation dataframe
  annotation_df <- data.frame(
    Gene = genes_df$names_plot,
    .y. = "Expression",
    group1 = control_group,
    group2 = test_group,
    y.position = genes_df$max + (minMaxNorm(genes_df$max, min(genes_df$max), max(genes_df$max))/2),
    p.adj = genes_df$padj,
    label = genes_df$signif
  )
  
  # Adjusting rownames
  rownames(annotation_df) <- annotation_df$Gene
  
  # Creating adjustment var for better plotting space
  countsMelted_df$adjustment <- countsMelted_df$Expression + (minMaxNorm(countsMelted_df$Expression, min(countsMelted_df$Expression), max(countsMelted_df$Expression)))
  
  # Ordering groups for plotting
  countsMelted_df[,group_var] <- factor(countsMelted_df[,group_var], levels = c(eval(control_group), eval(test_group)))
  
  # Setting gene order if available
  if (is.null(plot_order) == FALSE) {
    annotation_df$Gene <- factor(annotation_df$Gene, levels = plot_order)
    countsMelted_df$Gene <- factor(countsMelted_df$Gene, levels = plot_order)
  } else {
    message("Gene order was not provided. Proceeding with default order.")
  }
  
  # Plotting
  print(ggplot(data = countsMelted_df, aes_string(x = group_var, y="Expression", fill = group_var)) +
          geom_violin(trim = TRUE) +
          geom_point(aes(y = adjustment), alpha = 0) +
          facet_wrap( ~ Gene, scale = "free") +
          stat_summary(geom = "crossbar", fun = median, width = 0.4, 
                       size = 0.3, position = position_dodge(width=0.9), 
                       show.legend = FALSE) +
          scale_fill_manual(breaks = c(eval(control_group), eval(test_group)),
                            values = c(control_color, test_color)) +
          theme_bw()+
          theme(axis.text.x=element_blank(), 
                axis.title.x=element_blank(),
                axis.text.y = element_text(size = 12),
                axis.title.y = element_text(size = 14, face = "bold"),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.text.x = element_text(size = 14, face = "bold"),
                legend.position = "right",
                plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                plot.subtitle = element_text(hjust = 0.5, size = 10)) +
          labs(y = "Gene expression (log2)") +
          geom_bracket(inherit.aes = FALSE,
                       data = annotation_df, vjust = -0.1, size = 0.7, 
                       tip.length = 0, label.size = 4, 
                       aes(y.position = y.position, label = label, 
                           xmin = group1, xmax = group2)) -> violin_plot
  )
  
  return(list(data = genes_df, 
              plot_annotation = annotation_df, 
              violin_plot = violin_plot))
}