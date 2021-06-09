#  @param uv A data.frame with columns "id", "U", "deltaU", "V", and "deltaV".
#
#  @param ti A data.frame containing the baseline covariates. Note that
#                there is not an "id" column. Can be NULL.
#
#  @param td A data.frame containing the time dependent covariates. Note that
#    there is not an "id" column. Can be NULL.
#
#  @param uniqueEvent A numeric vector. The unique event times.
#
#  @param IPCW A numeric vector. The inverse probability weighting.
#
#' @importFrom survival coxph survfit
#' @importFrom stats as.formula
#' @include zooStep.R
.VStep <- function(uv, ti, td, uniqueEvent, IPCW) {

  n <- nrow(x = uv)

  # unique times to treatment discontinuation, clinical event, or censoring
  fit <- survival::coxph(formula = Surv(V, deltaV) ~ 1, data = uv)
  ss <- survival::survfit(formula = fit)

  uniqueV <- ss$time
  nUniqueV <- length(x = uniqueV)

  cumu1.hazard <- -log(x = ss$surv)

  # {nUniqueV}
  hazard.KM <- c(cumu1.hazard[1L], diff(x = cumu1.hazard))

  # {n x nUniqueV}
  hazard2.KM <- matrix(data = hazard.KM, 
                       nrow = n,  
                       ncol = nUniqueV,  
                       byrow = TRUE)

  # {n x nUniqueV}
  Kt.KM <- matrix(data = ss$surv, nrow = n, ncol = nUniqueV, byrow = TRUE)

  # {n x nUniqueV}
  ft.KM <- Kt.KM*hazard2.KM

  # {n x nUniqueV}
  vvEq <- outer(X = uv$V, Y = uniqueV-1e-8, FUN = ">") *
          outer(X = uv$V, Y = uniqueV+1e-8, FUN = "<")

  nVE <- nUniqueV + length(x = uniqueEvent)

  jjorder <- sort(x = c(uniqueV, uniqueEvent), index.return = TRUE)
  veTimes <- jjorder$x
  veTimesSortedID <- jjorder$ix

  # {n x nUniqueV}
  temp1 <- Kt.KM * hazard2.KM * vvEq
  temp2 <- rowSums(x = temp1)
  thetaVi <- matrix(data = temp2, nrow = n, ncol = nUniqueV)

  # {n x nVE}
  thetaVi.forD <- matrix(data = temp2, nrow = n, ncol = nVE)

  # {n x nUniqueV}
  temp1 <- Kt.KM * vvEq
  temp2 <- rowSums(x = temp1)
  KVi.KM <- matrix(data = temp2, nrow = n, ncol = nUniqueV)

  # {n x nUniqueV}
  thataBart <- Kt.KM
  # {n x nVE}
  thataBart.forD <- matrix(data = NA, nrow = n, ncol = nVE)
  thataBart.forD[,veTimesSortedID <= nUniqueV] <- thataBart
  thataBart.forD <- zooStep(mat = thataBart.forD)

  
  # {n*nUniqueV}
  VV <- rep(x = uv$V, each = nUniqueV)
  DV <- rep(x = uv$deltaV, each = nUniqueV)

  if (!is.null(x = ti)) {
    Lti <- apply(X = ti, MARGIN = 2L, FUN = rep, each = nUniqueV)
  } else {
    Lti <- NULL
  }

  start <- rep(x = c(0.0, uniqueV[-nUniqueV]), times = n)
  stop <- rep(x = uniqueV, times = n)


  dataTD <- data.frame("start" = start,
                       "stop" = stop,
                       "V.time" = {{VV > start} & {VV < stop+1e-8}} * DV)

  if (!is.null(x = td)) dataTD <- cbind(dataTD, td)
  if (!is.null(x = Lti)) dataTD <- cbind(dataTD, Lti)

  retain <- VV > start
  dataTD <- dataTD[retain,]
  mns <- colMeans(x = dataTD[,-c(1L:3L),drop=FALSE])
  dataTD <- scale(x = dataTD, center = c(0,0,0,mns), scale = FALSE)
  dataTD <- as.data.frame(x = dataTD)

  if (!is.null(x = td) && !is.null(x = ti)) {
    tdCols <- 4L:{4L+ncol(x = td)-1L}
    tiCols <- {4L+ncol(x = td)}:ncol(x = dataTD)

    Ltd <- scale(x = td, center = colMeans(x = td[retain,,drop=FALSE]), scale = FALSE)
    Lti <- scale(x = ti, center = colMeans(x = Lti[retain,,drop=FALSE]), scale = FALSE)

    form <- stats::as.formula(object = paste0("Surv(start,stop,V.time)~",
                                              paste(colnames(x = ti), collapse="+"), "+", 
                                              paste(colnames(x = td), collapse="+")))

  } else if (!is.null(x = td)) {
    tdCols <- 4L:{4L+ncol(x = td)-1L}
    tiCols <- NULL

    Ltd <- scale(x = td, center = colMeans(x = td[retain,,drop=FALSE]), scale = FALSE)
    Lti <- NULL

    form <- stats::as.formula(object = paste0("Surv(start,stop,V.time)~",
                                              paste(colnames(x = td), collapse="+")))

  } else if (!is.null(x = ti)) {
    tdCols <- NULL
    tiCols <- 4L:ncol(x = dataTD)

    Ltd <- NULL
    Lti <- scale(x = ti, center = colMeans(x = Lti[retain,,drop=FALSE]),
                 scale = FALSE)

    form <- stats::as.formula(object = paste0("Surv(start,stop,V.time)~",
                                              paste(colnames(x = ti), collapse="+")))
  }

  fit <- survival::coxph(formula = form, data = dataTD)

  fit_result <- list("form" = form, "fit" = fit)

  gammahat <- fit$coefficients

  ss <- survival::survfit(formula = fit)

  # {nUniqueV}
  cumu1.hazard <- -log(ss$surv)
  base.hazard1 <- c(cumu1.hazard[1L], diff(x = cumu1.hazard))
  hazard1 <- base.hazard1

  # {n x nUniqueV}
  if (!is.null(x = td) && !is.null(x = ti)) {
    licoef <- match(x = colnames(x = Lti), table = names(x = gammahat))
    ldcoef <- match(x = colnames(x = Ltd), table = names(x = gammahat))

    piTtilde <- matrix(data = hazard1, nrow = n, ncol = nUniqueV, byrow = TRUE)*
                exp(x = drop(x = Lti %*% gammahat[licoef]) + 
                        matrix(data = Ltd %*% gammahat[ldcoef], 
                               nrow = n, ncol = nUniqueV, byrow = TRUE))
  } else if (!is.null(x = ti)) {
    licoef <- match(x = colnames(x = Lti), table = names(x = gammahat))
    piTtilde <- matrix(data = hazard1, nrow = n, ncol = nUniqueV, byrow = TRUE)*
                exp(x = drop(x = Lti %*% gammahat[licoef]))
  } else if (!is.null(x = td)) {
    ldcoef <- match(x = colnames(x = Ltd), table = names(x = gammahat))
    piTtilde <- matrix(data = hazard1, nrow = n, ncol = nUniqueV, byrow = TRUE)*
                exp(x = matrix(data = Ltd %*% gammahat[ldcoef], 
                               nrow = n, ncol = nUniqueV, byrow = TRUE))
  }
  hazard2 <- piTtilde

  hazard2[hazard2 > 1.0] <- 1.0

  # {n x nUniqueV}
  temp1 <- 1.0 - piTtilde
  temp1[temp1 < 0.0] <- 1.0
  Kt <- t(x = apply(X = temp1, MARGIN = 1L, FUN = cumprod))

  # {n x nUniqueV}
  temp1 <- Kt*hazard2*vvEq
  temp2 <- rowSums(x = temp1)
  # {n x nVE}

  fVi.forD <- matrix(data = temp2, nrow = n, ncol = nVE)
  fVi.forD[fVi.forD <= 1e-8] <- 1.0

  Kt.forD <- matrix(data = NA, nrow = n, ncol = nVE)
  Kt.forD[,veTimesSortedID <= nUniqueV] <- Kt
  Kt.forD <- zooStep(mat = Kt.forD)

  w0 <- {outer(X = uv$V, Y = veTimes, FUN = ">")*uv$deltaV + {1L-uv$deltaV}} * 
        thataBart.forD / Kt.forD

  w1 <- {outer(X = uv$V, Y = veTimes+1e-8, FUN = "<")}*uv$deltaV*thetaVi.forD/fVi.forD

  w <- w0 + w1

  w <- w*IPCW

  # {n x nUniqueEvent}
  w <- w[,veTimesSortedID > nUniqueV]

  nUniqueEvent <- length(x = uniqueEvent)

  start <- rep(x = c(0.0, uniqueEvent[-nUniqueEvent]), times = n)
  stop <- rep(x = uniqueEvent, times = n)
  VV <- rep(x = uv$V, each = nUniqueEvent)
  UU <- rep(x = uv$U, each = nUniqueEvent)
  DD <- rep(x = uv$deltaU, each = nUniqueEvent)

  dataU <- data.frame("start" = start,
                      "stop" = stop,
                      "weight" = pmin(c(t(x = w)), 100.0),
                      "Z" = as.integer(x = VV > {stop - 1e-8}),
                      "U.time" = {{UU > start} & UU < {stop+1e-8}}*DD) 

  retain <- {UU > start} & {dataU$weight > 1e-8}
  dataU <- dataU[retain,]

  fit <- survival::coxph(formula = Surv(start, stop, U.time) ~ Z, 
                         weights = dataU$weight, 
                         data = dataU)

  beta <- fit$coefficients
  est <- exp(x = beta)
  est <- as.numeric(x = est)

  return( list("psi" = est, "coxPH" = fit_result) )

}

