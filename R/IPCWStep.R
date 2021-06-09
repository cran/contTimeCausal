#  @param uv A data.frame with columns "id", "U", "deltaU", "V", and "deltaV".
#
#  @param ti A data.frame containing the baseline covariates. Note that
#    there is not an "id" column. Can be NULL.
#
#' @importFrom survival coxph survfit
#' @importFrom stats as.formula
#' @include zooStep.R
.IPCWStep <- function(uv, ti) {

  # this gives us the unique times at which an event occurred

  fit <- survival::coxph(formula = Surv(U, deltaU) ~ 1, data = uv)
  ss <- survival::survfit(formula = fit)

  uniqueEvent <- ss$time[ss$n.event != 0L]
  nUniqueEvent <- length(x = uniqueEvent)

  # this gives us the unique times at which censoring occurred

  fit <- survival::coxph(formula = Surv(U, 1L-deltaU) ~ 1, data = uv)
  ss <- survival::survfit(formula = fit)

  uniqueCensor <- ss$time[ss$n.event != 0L]
  nUniqueCensor <- length(x = uniqueCensor)

  n <- nrow(x = uv)

  # time intervals repeated for each participant
  start <- rep(x = c(0.0, uniqueCensor[-nUniqueCensor]), times = n)
  stop <- rep(x = uniqueCensor, times = n)

  # if treatment discontinuation was optional, use V (U>V)
  # if treatment discontinuation was due to clinical event or censoring, 
  #   (U == V) shift V away from U
  VV <- rep(x = uv$V*uv$deltaV+(uv$V+1.0)*(1L-uv$deltaV), 
            each = nUniqueCensor)

  # participant's clinical event or censoring time repeated for each
  # time interval
  UU <- rep(x = uv$U, each = nUniqueCensor)

  # participant's indicator of clinical event (1) or censoring (0) repeated
  # for each time interval.
  DD <- rep(x = uv$deltaU, each = nUniqueCensor)

  if (!is.null(x = ti)) {
    # note this is now a matrix {n*nUnique x nX)
    # participant's baseline covariates repeated for each time point
    Lti <- apply(X = ti, MARGIN = 2L, FUN = rep, each = nUniqueCensor)
  } else {
    Lti <- NULL
  }

  # S.time is the indicator of censoring occurring during the time interval
  # Z is the indicator of treatment discontinuation, clinical event,
  # censoring, or treatment-terminating event occurring after the interval
  dataTS <- data.frame("start" = start,
                       "stop" = stop,
                       "S.time" = {{UU > start} & {UU < {stop+1e-8}}}*{1L-DD},
                       "Z" = {VV > stop}*1L)

  # if baseline covariates included in model, add to data.frame
  if (!is.null(x = Lti)) dataTS <- cbind(dataTS, Lti)

  # keep only those time intervals for which U > lower
  retain <- UU > start
  dataTS <- dataTS[retain,]

  if (!is.null(x = ti)) {

    # center the baseline covariates using the truncated data set
    mns <- colMeans(x = dataTS[ ,-(1L:4L) ,drop = FALSE])

    # center baseline covariates
    Lti <- scale(x = ti, center = mns, scale = FALSE)
    dataTS <- scale(x = dataTS, center = c(0.0,0.0,0.0,0.0,mns), scale = FALSE)

    dataTS <- as.data.frame(x = dataTS)

    # create formula for Cox regression
    form <- stats::as.formula(object = paste0("Surv(start, stop, S.time)~Z+",
                                              paste(colnames(x = Lti), collapse = "+")))
  } else {

    # create formula for Cox regression
    form <- stats::as.formula(object = paste0("Surv(start, stop, S.time)~Z"))

  }

  fit <- survival::coxph(formula = form, data = dataTS)
  gammahat <- fit$coefficients

  ss <- survival::survfit(formula = fit)

  # {nUniqueCensor}
  cumu3.hazard <- -log(x = ss$surv)
  hazard3 <- c(cumu3.hazard[1L], diff(x = cumu3.hazard))
  if (length(x = hazard3) > nUniqueCensor) {
    hazard3 <- c(0.0, hazard3)
  }


  if (!is.null(x = ti)) {
    # match column headers of Lti to estimated parameter names
    ltiCoef <- match(x = colnames(x = Lti), table = names(x = gammahat))

    # {n x nUniqueCensor}
    hazard3 <- matrix(data = hazard3, 
                      nrow = n, 
                      ncol = nUniqueCensor, 
                      byrow = TRUE) *
               exp(x = drop(Lti %*% gammahat[ltiCoef]) + 
                       outer(X = uv$V, 
                             Y = uniqueCensor-1e-8,  
                             FUN = ">")*gammahat[-ltiCoef])
  } else {

    hazard3 <- matrix(data = hazard3, 
                      nrow = n, 
                      ncol = nUniqueCensor, 
                      byrow = TRUE) *
               exp(x = outer(X = uv$V, 
                             Y = uniqueCensor-1e-8,  
                             FUN = ">")*gammahat)

  }

  # {n x nUniqueCensor}
  temp1 <- 1.0 - hazard3
  temp1[temp1 <= 1e-8] <- 1.0

  # {n x nUniqueCensor} not in ctsftm
  K3t <- t(x = apply(X = temp1, MARGIN = 1L, FUN = cumprod))

  # nUniqueCensor + nUniqueEvent === nCE
  nCE <- nUniqueCensor + nUniqueEvent

  jjorder <- sort(x = c(uniqueCensor,uniqueEvent), index.return = TRUE)
  ceTimesSortedID <- jjorder$ix

  # {n x nUniqueEvent}
  tmp <- outer(X = uv$U, Y = uniqueEvent-1e-8, FUN = ">") *
         outer(X = uv$U, Y = uniqueEvent+1e-8, FUN = "<")

  K3t.forD <- matrix(data = NA, nrow = n, ncol = nCE)
  K3t.forD[,ceTimesSortedID <= nUniqueCensor] <- K3t
  K3t.forD <- zooStep(mat = K3t.forD)
  K3t <- K3t.forD[,ceTimesSortedID > nUniqueCensor]

  # {n}
  K3 <- rowSums(x = tmp*K3t)*uv$deltaU

  IPCW <- numeric(length = n)
  IPCW[uv$deltaU == 1L] <- 1.0/K3[uv$deltaU == 1L]

  return( list("IPCW" = IPCW, "uniqueEvent" = uniqueEvent) )

}

