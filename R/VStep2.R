#  @param uv A data.frame with columns "id", "U", "deltaU", "V", and "deltaV".
#
#  @param ti A data.frame containing the baseline covariates. Note that
#                there is not an "id" column. Can be NULL.
#
#  @param td A data.frame containing the time dependent covariates. Note that
#    there is not an "id" column. Can be NULL.
#
#  @param uniqueEvent A numeric vector. The unique event times
#
#  @param IPCW A numeric vector. The IPCW weights.
#
#' @importFrom survival coxph survfit
#' @importFrom stats as.formula lm
.VStep2 <- function(uv, ti, td, uniqueEvent, IPCW) {

  n <- nrow(x = uv)

  # unique times to treatment discontinuation, clinical event, or censoring

  fit <- survival::coxph(formula = Surv(V, deltaV) ~ 1, data = uv)
  ss <- survival::survfit(formula = fit)

  uniqueV <- ss$time
  nUniqueV <- length(x = uniqueV)

  # repeat V and deltaV for each time interval
  VV <- rep(x = uv$V, each = nUniqueV)
  DV <- rep(x = uv$deltaV, each = nUniqueV)

  # repeat the time interval for each individual
  start <- rep(x = c(0.0, uniqueV[-nUniqueV]), times = n)
  stop <- rep(x = uniqueV, times = n)

  if (!is.null(x = ti)) {
    # repeat the baseline covariates for each time interval
    Lti <- apply(X = ti, MARGIN = 2L, FUN = rep, each = nUniqueV)
  } else {
    Lti <- NULL
  }

  dataTD <- data.frame("start" = start,
                       "stop" = stop,
                       "V.time" = {{VV > start} & {VV < stop+1e-8}} * DV)
  if (!is.null(x = td)) dataTD <- cbind(dataTD, td)
  if (!is.null(x = Lti)) dataTD <- cbind(dataTD, Lti)

  # retrain only those time intervals V > lower
  retain <- VV > start
  dataTD <- dataTD[retain,]

  mns <- colMeans(x = dataTD[,-c(1L:3L),drop=FALSE])

  # center covariates based on mean of truncated data
  dataTD <- scale(x = dataTD, center = c(0,0,0,mns), scale = FALSE)
  dataTD <- as.data.frame(x = dataTD)


  if (!is.null(x = td) & !is.null(x = ti)) {
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

  } else  if (!is.null(x = ti)) {
    tdCols <- NULL
    tiCols <- 4L:ncol(x = dataTD)

    Ltd <- NULL
    Lti <- scale(x = ti, center = colMeans(x = Lti[retain,,drop=FALSE]), scale = FALSE)

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
  if (!is.null(x = td) & !is.null(x = ti)) {
    ltiCoef <- match(x = colnames(x = Lti), table = names(x = gammahat))
    ldiCoef <- match(x = colnames(x = Ltd), table = names(x = gammahat))
    piTtilde <- matrix(data = hazard1, nrow = n, ncol = nUniqueV, byrow = TRUE)*
                exp(x = drop(x = Lti %*% gammahat[ltiCoef]) +  
                        matrix(data = Ltd %*% gammahat[ldiCoef], 
                               nrow = n, ncol = nUniqueV, byrow = TRUE))
    c4 <- exp(x = -(drop(x = Lti %*% gammahat[ltiCoef]) +  
                         matrix(data = Ltd %*% gammahat[ldiCoef], 
                                nrow = n, ncol = nUniqueV, byrow = TRUE)))
  } else if (!is.null(x = td)) {
    piTtilde <- matrix(data = hazard1, nrow = n, ncol = nUniqueV, byrow = TRUE)*
                exp(x = matrix(data = Ltd %*% gammahat, 
                               nrow = n, ncol = nUniqueV, byrow = TRUE))
    c4 <- exp(x = -matrix(data = Ltd %*% gammahat, 
                          nrow = n, ncol = nUniqueV, byrow = TRUE))
  } else if (!is.null(x = ti)) {
    piTtilde <- matrix(data = hazard1, nrow = n, ncol = nUniqueV, byrow = TRUE)*
                exp(x = drop(x = Lti %*% gammahat))
    c4 <- exp(x = -drop(x = Lti %*% gammahat))
  }

  hazard2 <- piTtilde

  hazard2[hazard2 > 1.0] <- 1.0

  # {n x nUniqueV}
  temp1 <- 1.0 - piTtilde
  temp1[temp1 < 0.0] <- 1.0
  Kt <- t(x = apply(X = temp1, MARGIN = 1L, FUN = cumprod))

  # I(V_i == t_j) Gamma
  dNu <- outer(X = uv$V, Y = uniqueV+1e-8, FUN = "<") * 
         outer(X = uv$V, Y = uniqueV-1e-8, FUN = ">") * uv$deltaV
  # I(V_i >= t_j)
  iXs1 <- outer(X = uv$V, Y = uniqueV-1e-8, FUN = ">")

  # I(U_i >= t_j)
  iXs2 <- outer(X = uv$U, Y = uniqueV-1e-8, FUN = ">")

  dMu <- dNu - piTtilde*iXs1*iXs2

  # IPCW is 0 when deltaU == 1 (clinical event occurred)
  use <- IPCW > 0.0

  if (!is.null(x = td) && !is.null(x = ti)) {

    tdn <- td[seq(from = 1L, to = n*nUniqueV, by = nUniqueV),,drop=FALSE]

    data1 <- data.frame("Y" = uv$U*{1.0-uv$deltaV} + uv$V*uv$deltaV,
                        tdn, ti)

    form <- as.formula(paste0("Y~", paste(colnames(x = ti), collapse = "+"),
                              "+", paste(colnames(x = td), collapse = "+")))

    lm.out <- stats::lm(formula = form, 
                        weights = IPCW[use], data = data1[use,])
    coeff <- lm.out$coefficients

    ltiCoef <- match(x = colnames(x = ti), table = names(x = coeff))
    ldiCoef <- match(x = colnames(x = tdn), table = names(x = coeff))

    Epart1 <- coeff[1L] + data.matrix(frame = ti) %*% coeff[ltiCoef] + 
                          data.matrix(frame = tdn) %*% coeff[ldiCoef]


    sum41c <- sum(c4*dMu*iXs1*iXs2*IPCW*drop(x = data1$Y-Epart1))

    data1$Y <- {{uv$U-uv$V}*uv$deltaV}
    lm.out <- stats::lm(formula = form, weights = IPCW[use], data = data1[use,])
    coeff <- lm.out$coefficients

    ltiCoef <- match(x = colnames(x = ti), table = names(x = coeff))
    ldiCoef <- match(x = colnames(x = tdn), table = names(x = coeff))

    Epart2 <- coeff[1L] + data.matrix(frame = ti) %*% coeff[ltiCoef] + 
                          data.matrix(frame = tdn) %*% coeff[ldiCoef]

    sum42c <- sum(c4*dMu*iXs1*iXs2*IPCW*drop(x = data1$Y-Epart2))

    est <- -sum42c/sum41c

  } else if (!is.null(x = td)) {

    tdn <- td[seq(from = 1L, to = n*nUniqueV, by = nUniqueV),,drop=FALSE]

    data1 <- data.frame("Y" = uv$U*{1.0-uv$deltaV} + uv$V*uv$deltaV,
                        tdn)

    form <- as.formula(paste0("Y~", paste(colnames(x = td), collapse = "+")))

    lm.out <- stats::lm(formula = form, 
                        weights = IPCW[use], data = data1[use,])
    coeff <- lm.out$coefficients

    ldiCoef <- match(x = colnames(x = tdn), table = names(x = coeff))

    Epart1 <- coeff[1L] + data.matrix(frame = tdn) %*% coeff[ldiCoef]

    sum41c <- sum(c4*dMu*iXs1*iXs2*IPCW*drop(x = data1$Y-Epart1))

    data1$Y <- {{uv$U-uv$V}*uv$deltaV}
    lm.out <- stats::lm(formula = form, weights = IPCW[use], data = data1[use,])
    coeff <- lm.out$coefficients

    ldiCoef <- match(x = colnames(x = tdn), table = names(x = coeff))

    Epart2 <- coeff[1L] + data.matrix(frame = tdn) %*% coeff[ldiCoef]

    sum42c <- sum(c4*dMu*iXs1*iXs2*IPCW*drop(x = data1$Y-Epart2))

    est <- -sum42c/sum41c


  } else if (!is.null(x = ti)) {

    data1 <- data.frame("Y" = uv$U*{1.0-uv$deltaV} + uv$V*uv$deltaV,
                        ti)

    form <- as.formula(paste0("Y~", paste(colnames(x = ti), collapse = "+")))

    lm.out <- stats::lm(formula = form, 
                        weights = IPCW[use], data = data1[use,])
    coeff <- lm.out$coefficients

    ltiCoef <- match(x = colnames(x = ti), table = names(x = coeff))

    Epart1 <- coeff[1L] + data.matrix(frame = ti) %*% coeff[ltiCoef]

    sum41c <- sum(c4*dMu*iXs1*iXs2*IPCW*drop(x = data1$Y-Epart1))

    data1$Y <- {{uv$U-uv$V}*uv$deltaV}
    lm.out <- stats::lm(formula = form, weights = IPCW[use], data = data1[use,])
    coeff <- lm.out$coefficients

    ltiCoef <- match(x = colnames(x = ti), table = names(x = coeff))

    Epart2 <- coeff[1L] + data.matrix(frame = ti) %*% coeff[ltiCoef]

    sum42c <- sum(c4*dMu*iXs1*iXs2*IPCW*drop(x = data1$Y-Epart2))

    est <- -sum42c/sum41c

  }

  return( list("psi" = est, "coxPH" = fit_result) )

}
