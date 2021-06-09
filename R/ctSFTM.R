#' Continuous-time Structural Failure Time Model
#'
#' The function estimates the regime (in terms of time to treatment initiation)
#'   of treatment effect for a survival outcome under a 
#'   Structural Failure Time Model (SFTM)
#'   with time-varying confounding in the presence of dependent censoring.
#'   Studying the effect of time to treatment discontinuation is applicable 
#'   by redefining "treatment initiation" in the current description to 
#'   "treatment discontinuation". 
#'
#' The SFTM assumes that the potential failure time \code{U} had the individual 
#'   never received treatment and the observed failure time \code{T} follow
#'   \deqn{U \sim \int_0^T e^{\psi A_u}d u, }
#'   where \code{~} means "has the same distribution as", and \eqn{A_u} is the 
#'   treatment indicator at time \eqn{u}.
#'   We assume that the individual continuously received treatment until 
#'   time \eqn{V}. The observed failure time can be censored assuming the 
#'   censoring time is independent of the failure time given the treatment and 
#'   covariate history (the so-called ignorable censoring). The current 
#'   function allows for multi-dimensional baseline covariates and/or
#'   multi-dimensional time-dependent covariate.
#'   Variance estimates should be implemented by delete-one-group jackknifing 
#'   and recalling ctSFTM.
#'
#'   If only time-independent covariates are included, the data.frame must 
#'   contain the following columns:
#'   \describe{
#'     \item{id:}{A unique participant identifier.}
#'     \item{U:}{The time to the clinical event or censoring.}
#'     \item{deltaU:}{The clinical event indicator (1 if U is the event time;
#'                    0 otherwise.}
#'     \item{V:}{The time to optional treatment discontinuation, a clinical 
#'               event, censoring, or a treatment-terminating event.}
#'     \item{deltaV:}{The indicator of optional treatment discontinuation 
#'                    (1 if treatment discontinuation was optional; 0 if 
#'                    treatment discontinuation was due to a clinical event, 
#'                    censoring or a treatment-terminating event.}
#'   }
#'
#' If time-dependent covariates are to be included, the data.frame must be
#'   a time-dependent dataset as described by package survival. Specifically,
#'   the time-dependent data must be specified for an interval (lower,upper]
#'   and the data must include the following additional columns:
#'   \describe{
#'      \item{start:}{The lower boundary of the time interval to which the
#'                    data pertain.}
#'      \item{stop:}{The upper boundary of the time interval to which the
#'                    data pertain.}
#'   }
#'
#'
#' @param data A data.frame object. A data.frame containing all observed data.
#'    At a minimum, this data.frame must contain columns with headers 
#'    "id", "U", "V", "deltaU", and "deltaV". If time-dependent covariates are
#'    included, additional columns include "stop" and "start". See Details for
#'    further information
#'
#' @param base A character or integer vector or NULL. The columns of data to be
#'   included in the time-independent component of the model. If NULL, 
#'   time-independent covariates are excluded from the Cox model for 
#'   treatment discontinuation. 
#'
#' @param td A character or integer vector or NULL. The columns of data to be
#'   included in the time-dependent component of the model. If NULL, 
#'   time-dependent covariates are excluded from the Cox model for 
#'   treatment discontinuation. 
#'
#' @returns An S3 object of class ctc. Object contains element `psi', the 
#'   estimate of the SFTM parameter(s) and `coxPH', the Cox 
#'   regression for V.
#'
#' @seealso \code{\link{ctCoxMSM}}
#'
#' @references 
#'   Yang, S., K. Pieper, and F. Cools. (2019) 
#'   Semiparametric estimation of structural failure time model in 
#'   continuous-time processes.
#'   Biometrika, 107(1), 123-136.
#'
#' @examples
#'
#'   data(ctcData)
#'
#'  # sample data to reduce computation time of example
#'  smp <- ctcData$id %in% sample(1:1000, 200, FALSE)
#'  ctcData <- ctcData[smp,]
#'
#'  # analysis with both time-dependent and time-independent components
#'  res <- ctSFTM(data = ctcData, base = "x", td = "xt")
#'
#'  # analysis with only the time-independent component
#'  res <- ctSFTM(data = ctcData, base = "x")
#'
#'  # analysis with only the time-dependent component
#'  res <- ctSFTM(data = ctcData, td = "xt")
#'
#' @export
#'
#' @include verifyInputs.R IPCWStep.R VStep2.R
ctSFTM <- function(data, base = NULL, td = NULL) {

  inputs <- .verifyInputs(df = data, base = base, td = td)

  ss <- .IPCWStep(uv = inputs$uv,
                  ti = inputs$ti)

  vstep <- .VStep2(uv = inputs$uv,
                   ti = inputs$ti,  
                   td = inputs$td,
                   uniqueEvent = ss$uniqueEvent, 
                   IPCW = ss$IPCW)

  message("Estimated psi: ", round(x = vstep$psi, digits = 4L))

  class(vstep) <- "ctc"

  return( vstep )

}
