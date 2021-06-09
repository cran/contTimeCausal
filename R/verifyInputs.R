# @param df A data.frame object. The data for the analysis. Can be a time 
#   dependent data set. It is assumed that this data.frame
#   contains column headers "U", "V", "deltaU", "deltaV", "id". If time
#   dependent covariates are included additional required columns are "start"
#   and "stop". These column headers should not be included in the base or td
#   inputs, but a check is in place to remove them if they are.
#
# @param base A character or numeric vector indicating the columns
#   of df corresponding to the time independent covariates. Can by NULL
#   indicating that no baseline covariates are included.
#
# @param td A character or numeric vector indicating the columns of df
#   corresponding to the time dependent covariates. Can by NULL
#   indicating that no time dependent covariates are included.
#
# @returns a list containing
#   \describe{
#      \item{uv}{A data.frame with columns "id", "U", "deltaU", "V", and 
#                "deltaV".}
#      \item{ti}{A data.frame containing the baseline covariates. Note that
#                there is not an "id" column.}
#      \item{td}{A data.frame containing the time dependent covariates
#                rebinned according to the unique V values. Note that
#                there is not an "id" column.}
#   }
#
#' @import dplyr
.verifyInputs <- function(df, base, td) {

  if (!is.data.frame(x = df)) {
    df <- as.data.frame(x = df)
  }

  if (any(is.na(x = df))) {
    stop("data must be complete", call. = FALSE)
  }

  # must provide at least one of inputs base and td
  if (is.null(x = base) && is.null(x = td)) {
    stop("at least one of inputs base and td must be specified", call. = FALSE)
  }

  # columns that make up the "uv" internal data.frame
  required <- c("id", "V", "deltaV", "U", "deltaU")

  # data.frame must contain columns with headers
  # U, V, deltaU, deltaV, id
  if (!all(required %in% colnames(x = df))) {
    stop("data.frame must contain column headers ", paste(required, collapse=", "),
         call. = FALSE)
  }

  # create "uv" data.frame and remove duplicate rows
  df_uv <- df[,required]
  df_uv <- df_uv %>% dplyr::distinct()

  # ensure that deltaV is integer or can be converted to integer without
  # loss of information
  if (!is.integer(x = df_uv$deltaV)) {
    if (is.logical(x = df_uv$deltaV)) {
      df_uv$deltaV <- as.integer(x = df_uv$deltaV)
    } else {
      tmp <- as.integer(x = round(x = df_uv$deltaV, digits = 0L))
      if (!isTRUE(x = all.equal(target = tmp, current = df_uv$deltaV))) {
        stop("deltaV must be integer", call. = FALSE)
      }
      df_uv$deltaV <- tmp
    }
  }

  # ensure that deltaV is binary 0/1
  if (!all(df_uv$deltaV %in% c(0L,1L))) {
    stop("deltaV must be integer 0/1", call. = FALSE)
  }

  # ensure that deltaU is integer or can be converted to integer without
  # loss of information
  if (!is.integer(x = df_uv$deltaU)) {
    if (is.logical(x = df_uv$deltaU)) {
      df_uv$deltaU <- as.integer(x = df_uv$deltaU)
    } else {
      tmp <- as.integer(x = round(x = df_uv$deltaU, digits = 0L))
      if (!isTRUE(x = all.equal(target = tmp, current = df_uv$deltaU))) {
        stop("deltaU must be integer", call. = FALSE)
      }
      df_uv$deltaU <- tmp
    }
  }

  # ensure that deltaU is binary 0/1
  if (!all(df_uv$deltaU %in% c(0L,1L))) {
    stop("deltaU must be integer 0/1", call. = FALSE)
  }


  # baseline covariates

  if (!is.null(x = base)) {
    # if base is a numeric/integer object, pull appropriate column headers
    if (!is.character(x = base)) {
      if (is.numeric(x = base)) {

        # ensure that all elements of base are > 0 and < # of columns
        if (max(base) > ncol(x = df) || any(base <= 0L)) {
          stop("inappropriate column index provided for input base", 
               call. = FALSE)
        }
        base <- colnames(x = df)[base]
      } else {
        stop("if provided, base must be a character or integer vector", 
             call. = FALSE)
      }
    }

    # ensure that none of the required headers are provided in base
    if (any(required %in% base)) {
      base <- base[!(base %in% required)]
      if (length(x = base) == 0L) {
        stop("inappropriate column index provided for input base", 
             call. = FALSE)
      }
    }

    # ensure that the provided column headers are found in the data
    if (!all(base %in% colnames(x = df))) {
      stop("baseline covariates ", 
           paste(base[!(base %in% colnames(x = df))], collapse = ", "),
           " not found in provided data.frame", call. = FALSE)
    }

    # extract baseline terms and remove duplicate rows and id column
    df_ti <- df[,c("id",base)]
    df_ti <- df_ti %>% dplyr::distinct()
    df_ti <- df_ti[, -1L, drop = FALSE]

    if (nrow(x = df_ti) != nrow(x = df_uv)) {
      stop("number of rows for time-independent covariates does not agree ",
           "with number of participants", call. = FALSE)
    }

  } else {
    df_ti <- NULL
  }

  #  time dependent covariates

  if (!is.null(x = td)) {

    # if td is a numeric/integer object, pull appropriate column headers
    if (!is.character(x = td)) {
      if (is.numeric(x = td)) {

        # ensure that all elements of td are > 0 and < # of columns
        if (max(td) > ncol(x = df) || any(td <= 0L)) {
          stop("inappropriate column index provided for input td", 
               call. = FALSE)
        }
        td <- colnames(x = df)[td]
      } else {
        stop("if provided, td must be a character or integer vector", 
             call. = FALSE)
      }
    }

    # include 'stop' and 'start' in required column headers
    required <- c(required, "stop", "start")

    # data.frame must contain columns with headers
    # U, V, deltaU, deltaV, id, start, stop
    if (!all(required %in% colnames(x = df))) {
      stop("data.frame must contain column headers ", required,
           call. = FALSE)
    }

    # ensure that none of the required headers are provided in td
    if (any(required %in% td)) {
      td <- td[!(td %in% required)]
      if (length(x = td) == 0L) {
        stop("inappropriate column index provided for input td", 
             call. = FALSE)
      }
    }

    # ensure that the provided column headers are found in the data
    if (!all(td %in% colnames(x = df))) {
      stop("time dependent covariates ", td[!(td %in% colnames(x = df))],
           " not found in provided data.frame", call. = FALSE)
    }

    # extract time dependent terms
    df_td <- df[,c("id", "start", "stop", td)]

    # retrieve the unique V times

    fit <- survival::coxph(Surv(V, deltaV) ~ 1 , data = df_uv)

    # the boundaries of each time interval
    upper <- survival::survfit(formula = fit)$time
    nt <- length(x = upper)
    lower <- c(0.0, upper[-nt])

    # matrix to hold the time dependent covariates for use in Vstep
    # {n*nt x ntd}
    tdOnV <- NULL

    for (i in 1L:nrow(x = df_uv)) {

      # identify the rows of time dependent data for participant i
      idCol <- which(x = df_td$id == df_uv$id[i])
      nTD <- length(x = idCol)

      cnt <- integer(length = nt) + length(x = idCol) + 1L

      for (j in 1L:nTD) {
        cnt <- cnt - {upper < {df_td$stop[idCol[j]] + 1e-8}}
      }

      vals <- data.matrix(frame = rbind(df_td[idCol,], df_td[idCol[nTD],]))
      tmp <- unname(obj = vals[cnt,][,td,drop=FALSE])

      tdOnV <- rbind(tdOnV, tmp)

    }

    tdOnV <- t(x = tdOnV)

    rownames(x = tdOnV) <- td
    td <- as.data.frame(x = t(x = tdOnV))
  } else {
    td <- NULL
  }

  return( list("uv" = df_uv, "ti" = df_ti, "td" = td) )

}
