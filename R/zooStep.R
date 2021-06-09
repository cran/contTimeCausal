# @param m x n matrix
#' @importFrom zoo zoo index na.approx
zooStep <- function(mat) {
  # { n x m+1}
  eg1 <- cbind(1L:ncol(x = mat), t(x = mat))

  # order the matrix
  egz1 <- zoo::zoo(x = eg1)

  # reset the index attribute to be first column of identifiers
  zoo::index(x = egz1) <- egz1[,1L]

  # replace NAs with interpolated results where the outside the interval the 
  # closest data extreme is used (rule = 2)
  temp1 <- zoo::na.approx(object = egz1, rule = 2L)

  # remove identifer column
  temp1 <- temp1[,-1L]

  # convert back to a matrix
  temp1 <- as.matrix(x = temp1)

  # return in the original dimension
  return( t(x = temp1) )

}

