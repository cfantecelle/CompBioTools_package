#' Check numeric(0)
#'
#' Checks if something returns numeric(0). Posted by 'Andrii' user at https://stackoverflow.com/questions/10664662/how-to-test-when-condition-returns-numeric0-in-r
#' @param x The value/array to check.
#' @return True or False.
#' @export
is.numeric0 <- function(x) {
  identical(x, numeric(0))
}

#' Min/max normalization
#'
#' Returns min/max normalized values.
#' @param x The vector of values to normalize according to the minimum and maximum values of the same vector.
#' @param min Minimum value in set of values. If min = NULL then min(x) is applied.
#' @param max Maximum value in set of values. If max = NULL then max(x) is applied.
#' @return Normalized values in the vector.
#' @export
minMaxNorm <- function(x, min = NULL, max = NULL) {
  
  if (is.null(min)) {
    min = min(x)
  }
  
  if (is.null(max)) {
    min = max(x)
  }
  
  y <- (x-min/(max-min))
        
  return(y)
  
}

#' DF to Matrix with rownames.
#'
#' Converts a dataframe to a matrix and keeps the rownames.
#' @param x The dataframe to be converted.
#' @return A matrix with the same rownames as the input dataframe.
#' @export
matrixWithRowNames <- function(x) {
  
  rows <- rownames(x)
  m <- as.matrix(x)
  rownames(m) <- rows
  return(m)
  
}

#' Get lower triangle
#'
#' Extracts the 'lower triangle' of a matrix.
#' @param x The original matrix.
#' @return The lower triangle of the original matrix.
#' @export
get_lower_tri<-function(x) {
  
  x[upper.tri(x)] <- NA
  
  return(x)
}

#' Get upper triangle
#'
#' Extracts the 'upper triangle' of a matrix.
#' @param x The original matrix.
#' @return The upper triangle of the original matrix.
#' @export
get_upper_tri <- function(x) {
  
  x[lower.tri(x)]<- NA
  
  return(x)
}

#' Not in operator
#'
#' The reverse of the `%in%` in operator. Tells if a value is not present in another set o values. 
#' @return True/False.
#' @export
`%notin%` <- Negate(`%in%`)