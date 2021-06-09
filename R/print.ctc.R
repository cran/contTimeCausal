#' Print the Primary Results
#'
#' Print the primary results of a ctCoxMSM() or ctSFTM() analysis.
#'
#' @param x An S3 object of class ctc. The value object returned by a call to 
#'   ctCoxMSM() or ctSFTM().
#'
#' @param ... ignored
#'
#' @export
#' @name print
#' @method print ctc
#'
#' @returns No return value, called to display key results.
#'
#' @examples
#'
#'  data(ctcData)
#'
#'  # sample data to reduce computation time of example
#'  smp <- ctcData$id %in% sample(1:1000, 150, FALSE)
#'  ctcData <- ctcData[smp,]
#'
#'  # analysis with both time-dependent and time-independent components
#'  res <- ctCoxMSM(data = ctcData, base = "x", td = "xt")
#'
#'  print(x = res)
#' 
#'  # analysis with both time-dependent and time-independent components
#'  res <- ctSFTM(data = ctcData, base = "x", td = "xt")
#'
#'  print(x = res)
#'
#' @importFrom methods show
print.ctc <- function(x, ...) {

  cf <- as.character(x = x$coxPH$form)
  message("form: ", cf[2L], " ~ ", cf[3L])
  show(x$coxPH$fit)
  
  cat("\n\nEstimated psi: ", round(x = x$psi, digits = 4L), "\n")

}

