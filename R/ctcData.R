#' Toy Dataset For Illustration
#'
#' This data set is provided for the purposes of illustrating the use of
#' the software. It includes a one-dimensional baseline covariate and a
#' one-dimensional time-dependent covariate.
#' 
#' 
#' @usage data(ctcData)
#'
#' @format ctcData is a data.frame containing data for 1,000 participants. The
#'   data.frame contains 9 columns: 
#'   \describe{
#'     \item{id}{An integer participant identifier.}
#'     \item{start}{The left side of the time interval for time-dependent covariate xt.}
#'     \item{stop}{The right side of the time interval for time-dependent covariate xt.}
#'     \item{xt}{A continuous time-dependent covariate.}
#'     \item{x}{A continuous baseline covariate.}
#'     \item{deltaU}{A binary indicator of the clinical event. If the
#'                   clinical event occurred, takes value 1; otherwise 0.}
#'     \item{deltaV}{A binary indicator of treatment discontinuation. If
#'                   treatment discontinuation was optional, takes value 1.
#'                   If treatment discontinuation was due to the clinical
#'                   event, censoring, or a treatment-terminating event, takes
#'                   value 0.}
#'     \item{U}{The time to the clinical event or censoring.}
#'     \item{V}{The time to optimal treatment discontinuation, the clinical
#'              event, censoring, or a treatment-terminating event.}
#'   }
#'
#' @name ctcData
#' @keywords datasets
NULL
