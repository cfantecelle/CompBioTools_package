#' Convert p-values to * codes.
#'
#' Converts p-values to standard * codes.
#' 
#' p < 0.05 = *
#' p < 0.01 = **
#' p < 0.001 = ***
#' p < 0.0001 = ****
#' 
#' @param annotation_df Annotation dataframe constaining stats values.
#' @param pvalue_col Name of the column that contains the p-values.
#' @return The same annotation_df with an extra column called 'pcode' with coded p-values.
#' @export
pvalueCodes <- function(annotation_df, pvalue_col) {
  
  annotation_df$pcode <- ""
  
  for(i in 1:length(annotation_df[[pvalue_col]])) {
    if (is.na(annotation_df[i,pvalue_col])) {
      annotation_df[i,"pcode"] <- NA
    } else if (annotation_df[i,pvalue_col] < 0.05 & annotation_df[i,pvalue_col] >= 0.01) {
      annotation_df[i,"pcode"] <- "*"
    } else if (annotation_df[i,pvalue_col] < 0.01 & annotation_df[i,pvalue_col] >= 0.001) {
      annotation_df[i,"pcode"] <- "**"
    } else if (annotation_df[i,pvalue_col] < 0.001 & annotation_df[i,pvalue_col] >= 0.0001){
      annotation_df[i,"pcode"] <- "***"
    } else if (annotation_df[i,pvalue_col] < 0.0001){
      annotation_df[i,"pcode"] <- "****"
    } else {
      annotation_df[i,"pcode"] <- ""
    }
  }
  
  return(annotation_df)
}

