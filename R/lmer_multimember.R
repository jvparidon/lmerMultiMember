#' @title multimembership linear mixed effects models
#' @description lme4::lmer but with multimembership random effects
#' @name lmer
#' @inheritParams lme4::lmer
#' @param memberships named list of weight matrices that will replace any
#' (dummy) random effects with matching names
#' @return lme4 model object
#' @export
#' @import lme4
#' @examples
#'
#' df <- data.frame(
#'   x = seq(60) + runif(60, 0, 10),
#'   y = seq(60) + rep(runif(6, 0, 10), 10),
#'   memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
#' )
#' weights <- weights_from_vector(df$memberships)
#'
#' # note that the grouping variable name is arbitrary -- it just has
#' # to match the name in the list and doesn't need to correspond to a column
#' # name in the data
#' lmer(y ~ x + (1 | members),
#'   data = df,
#'   memberships = list(members = weights)
#' )
lmer <- function(formula,
                 data = NULL,
                 REML = TRUE,
                 control = lme4::lmerControl(),
                 start = NULL,
                 verbose = 0L,
                 weights = NULL,
                 na.action = na.omit,
                 offset = NULL,
                 contrasts = NULL,
                 devFunOnly = FALSE,
                 memberships = NULL) {
  orig_call <- match.call()

  # detect if lmerTest is loaded
  if ("lmerTest" %in% (.packages())) {
    lmerlib <- "lmerTest"
  } else {
    lmerlib <- "lme4"
  }

  # if there are no multiple membership matrices in the function call
  # just pass it to lme4 or lmerTest instead
  if (is.null(memberships)) {
    return(do.call(get("lmer", asNamespace(lmerlib)), args = list(
      formula,
       data = data,
       REML = REML,
       control = control,
       start = start,
       verbose = verbose,
       weights = weights,
       na.action = na.action,
       offset = offset,
       contrasts = contrasts,
       devFunOnly = devFunOnly
    )))
  }

  # get names of multimembership variables
  multi_RE_names <- names(memberships)

  # subset data to include only variables used in formula
  # this is used for correct missing data handling
  data <- data[setdiff(all.vars(formula), multi_RE_names)]

  # get index of bars (location of random effects) in model formula
  bar_idx <- lme4::findbars(formula)

  # get random effect grouping variables
  RE_vars <- vapply(bar_idx, function(x) deparse(x[[3]]), character(1))

  # initialize list that will hold multimembership random effects matrices
  Ztlist <- list()

  # iterate over random effects
  for (i in seq_along(bar_idx)) {
    RE_name <- deparse(bar_idx[[i]])

    # find corresponding random-effects term
    RE_idx <- which(multi_RE_names == RE_vars[i])

    if (length(RE_idx) > 0) {
      # select relevant weight matrix
      M <- memberships[[RE_idx]][, complete.cases(data)]

      # extract LHS (effect)
      subformula <- as.formula(substitute(~z, list(z = bar_idx[[i]][[2]])))

      # construct model matrix & compute Khatri-Rao product
      X <- model.matrix(subformula, data[complete.cases(data), , drop=FALSE])
      Zt <- Matrix::KhatriRao(M, t(X), make.dimnames = TRUE)

      # FIXME: mess with names?
      Ztlist[[RE_name]] <- Zt

      # if necessary, add factor to data
      if (!RE_vars[i] %in% names(data)) {
        # if the factor has non-trivial ordering, it should be included
        # TODO: do we have to worry about ordering of Z? test!
        data[[RE_vars[i]]] <- rep_len(
          factor(rownames(memberships[[RE_idx]])),
          nrow(data)
        )
      }
    }
  }

  # create model specification but don't fit
  lmod <- lFormula(formula,
    data = data,
    REML = REML,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    control = control
  )

  # substitute new Ztlist elements into the model specification
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)

  # finish fitting
  devfun <- do.call(mkLmerDevfun, lmod)
  if (devFunOnly) {
    return(devfun)
  }

  opt <- optimizeLmer(devfun,
    optimizer = control$optimizer,
    restart_edge = control$restart_edge,
    boundary.tol = control$boundary.tol,
    start = NULL,
    verbose = 0L,
    control = control$optCtrl
  )

  res <- mkMerMod(
    environment(devfun),
    opt,
    lmod$reTrms,
    fr = lmod$fr,
    mc = orig_call
  )

  # convert model object to lmerTest or lmer, depending on what is loaded
  if (lmerlib == "lmerTest") {
    # first, convert model object to lmerModLmerTest
    res <- lmerTest:::as_lmerModLT(res, devfun)

    # convert model object to lmerModLmerTestMultiMember and add memberships
    res <- as(res, "lmerModLmerTestMultiMember")
    res@memberships <- memberships
  } else {
    # convert model object to lmerModMultiMember and add memberships
    res <- as(res, "lmerModMultiMember")
    res@memberships <- memberships
  }

  return(res)

}


#' @title multimembership generalized linear mixed effects models
#' @description lme4::glmer but with multimembership random effects
#' @name glmer
#' @inheritParams lme4::glmer
#' @param memberships named list of weight matrices that will replace any
#' (dummy) random effects with matching names
#' @return lme4 model object
#' @export
#' @import lme4
#' @examples
#'
#' df <- data.frame(
#'   x = runif(60, 0, 1),
#'   y = rbinom(60, 1, 0.6),
#'   memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
#' )
#' weights <- weights_from_vector(df$memberships)
#'
#' # note that the grouping variable name is arbitrary -- it just has
#' # to match the name in the list and doesn't need to correspond to a column
#' # name in the data
#' glmer(y ~ x + (1 | members),
#'   data = df,
#'   family = binomial,
#'   memberships = list(members = weights)
#' )
glmer <- function(formula,
                  data = NULL,
                  family,
                  control = lme4::glmerControl(),
                  start = NULL,
                  verbose = 0L,
                  nAGQ = 1L,
                  weights = NULL,
                  na.action = na.omit,
                  offset = NULL,
                  contrasts = NULL,
                  devFunOnly = FALSE,
                  memberships = NULL) {
  orig_call <- match.call()

  +  # if there are no multiple membership matrices in the function call
    +  # just pass it to lme4 instead
  if (is.null(memberships)) {
    return(lme4::glmer(formula,
      data = data,
      family = family,
      control = control,
      start = start,
      verbose = verbose,
      nAGQ = nAGQ,
      weights = weights,
      na.action = na.action,
      offset = offset,
      contrasts = contrasts,
      devFunOnly = devFunOnly
    ))
  }

  # get names of multimembership variables
  multi_RE_names <- names(memberships)

  # subset data to include only variables used in formula
  # this is used for correct missing data handling
  data <- data[setdiff(all.vars(formula), multi_RE_names)]

  # get index of bars (location of random effects) in model formula
  bar_idx <- lme4::findbars(formula)

  # get random effect grouping variables
  RE_vars <- vapply(bar_idx, function(x) deparse(x[[3]]), character(1))

  # initialize list that will hold multimembership random effects matrices
  Ztlist <- list()

  # iterate over random effects
  for (i in seq_along(bar_idx)) {
    RE_name <- deparse(bar_idx[[i]])

    # find corresponding random-effects term
    RE_idx <- which(multi_RE_names == RE_vars[i])

    if (length(RE_idx) > 0) {
      # select relevant weight matrix
      M <- memberships[[RE_idx]][, complete.cases(data)]

      # extract LHS (effect)
      subformula <- as.formula(substitute(~z, list(z = bar_idx[[i]][[2]])))

      # construct model matrix & compute Khatri-Rao product
      X <- model.matrix(subformula, data[complete.cases(data), , drop=FALSE])
      Zt <- Matrix::KhatriRao(M, t(X), make.dimnames = TRUE)

      # FIXME: mess with names?
      Ztlist[[RE_name]] <- Zt

      # if necessary, add factor to data
      if (!RE_vars[i] %in% names(data)) {
        # if the factor has non-trivial ordering, it should be included
        # TODO: do we have to worry about ordering of Z? test!
        data[[RE_vars[i]]] <- rep_len(
          factor(rownames(memberships[[RE_idx]])),
          nrow(data)
        )
      }
    }
  }

  # create model specification but don't fit
  glmod <- glFormula(formula,
    data = data,
    family = family,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    control = control
  )

  # substitute new Ztlist elements into the model specification
  for (m in names(Ztlist)) {
    glmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  glmod$reTrms$Zt <- do.call(rbind, glmod$reTrms$Ztlist)

  # finish fitting
  ##############################################################################
  ## TODO: figure out how much of the nAGQ wrangling can be handed off to lme4
  ## because the tangle of if/else statements below is pretty ugly
  ##############################################################################
  nAGQinit <- if (control$nAGQ0initStep) 0L else 1L
  devfun <- do.call(mkGlmerDevfun, c(glmod, list(
    verbose = verbose,
    control = control,
    nAGQ = nAGQinit
  )))
  if (nAGQ == 0 && devFunOnly) {
    return(devfun)
  }

  if (is.list(start)) {
    start.bad <- setdiff(names(start), c("theta", "fixef"))
    if (length(start.bad) > 0) {
      stop(sprintf(
        "bad name(s) for start vector (%s); should be %s and/or %s",
        paste(start.bad, collapse = ", "),
        shQuote("theta"),
        shQuote("fixef")
      ), call. = FALSE)
    }
    if (!is.null(start$fixef) && nAGQ == 0) {
      stop("should not specify both start$fixef and nAGQ==0")
    }
  }

  # FIXME: allow calc.derivs, use.last.params etc. if nAGQ=0
  if (control$nAGQ0initStep) {
    opt <- optimizeGlmer(devfun,
      optimizer = control$optimizer[[1]],
      # DON'T try fancy edge tricks unless nAGQ=0 explicitly set
      restart_edge = if (nAGQ == 0) control$restart_edge else FALSE,
      boundary.tol = if (nAGQ == 0) control$boundary.tol else 0,
      control = control$optCtrl,
      start = start,
      nAGQ = 0,
      verbose = verbose,
      calc.derivs = FALSE
    )
  }

  if (nAGQ > 0L) {

    # update deviance function to include fixed effects as inputs
    devfun <- updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)

    # updateStart is a utility function from lme4 that is necessary here
    # but for some reason lme4 doesn't export it
    updateStart <- function(start, theta) {
      if (is.numeric(start)) {
        theta
      } else {
        if (!is.null(start$theta)) {
          start$theta <- theta
        }
        start
      }
    }

    if (control$nAGQ0initStep) {
      start <- updateStart(start, theta = opt$par)
    }
    ## if nAGQ0 was skipped
    ## we don't actually need to do anything here, it seems --
    ## getStart gets called again in optimizeGlmer

    if (devFunOnly) {
      return(devfun)
    }
    # reoptimize deviance function over covariance parameters and fixed effects
    opt <- optimizeGlmer(devfun,
      optimizer = control$optimizer[[2]],
      restart_edge = control$restart_edge,
      boundary.tol = control$boundary.tol,
      control = control$optCtrl,
      start = start,
      nAGQ = nAGQ,
      verbose = verbose,
      stage = 2,
      calc.derivs = control$calc.derivs,
      use.last.params = control$use.last.params
    )
  }
  ##############################################################################

  m1 <- mkMerMod(environment(devfun),
    opt,
    glmod$reTrms,
    fr = glmod$fr,
    mc = orig_call
  )

  # convert model object to glmerModMultiMember and add memberships
  m1 <- as(m1, "glmerModMultiMember")
  m1@memberships <- memberships
  return(m1)
}


#' @importFrom lme4 fixef
#' @export
lme4::fixef

#' @importFrom lme4 ranef
#' @export
lme4::ranef


#' @title Model object for multimembership linear mixed models
#' @description The \code{merModMultiMember} class extends \code{merMod}
#' from the \pkg{lme4}-package.
#' @seealso \code{\link[lme4]{lmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @importClassesFrom lme4 merMod
#' @return An object of class \code{merModMultiMember} similar to
#' \code{merMod} objects (see \code{\link[lme4]{merMod}}) with
#' extra information about the multimembership random effects.
merModMultiMember <-
  setClass("merModMultiMember",
    contains = c("merMod"), slots = c(memberships = "list")
  )


#' @title Model object for multimembership linear mixed models
#' @description The \code{lmerModMultiMember} class extends \code{lmerMod}
#' from the \pkg{lme4}-package and \code{merModMultiMember}.
#' @seealso \code{\link[lme4]{lmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @importClassesFrom lme4 lmerMod
#' @return An object of class \code{lmerModMultiMember} similar to
#' \code{merModMultiMember} objects but inheriting from \code{lmerMod}
lmerModMultiMember <-
  setClass("lmerModMultiMember",
           contains = c("lmerMod", "merModMultiMember")
  )


#' @title Model object for multimembership linear mixed models with lmerTest
#' @description The \code{lmerModMultiMember} class extends
#' \code{lmerModLmerTest} from the \pkg{lmerTest}-package and
#' \code{lmerModMultiMember}.
#' @seealso \code{\link[lmerTest]{lmer}} and \code{\link[lme4]{lmer}}
#' @export
#' @return An object of class \code{lmerModMultiMember} similar to
#' \code{merModMultiMember} objects but inheriting from \code{lmerMod}
lmerModLmerTestMultiMember <-
  setClass(
    "lmerModLmerTestMultiMember",
    contains = c("lmerModMultiMember",
    #             `if`(("lmerTest" %in% (.packages())),
    #                    "lmerModLmerTest", NULL)
    `if`((require(lmerTest)), "lmerModLmerTest", NULL)
    )
  )


#' @title Model object for multimembership generalized linear mixed models
#' @description The \code{glmerModMultiMember} class extends \code{glmerMod}
#' from the \pkg{lme4}-package and \code{merModMultiMember}.
#' @seealso \code{\link[lme4]{glmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @importClassesFrom lme4 glmerMod
#' @return An object of class \code{glmerModMultiMember} similar to
#' \code{merModMultiMember} objects but inheriting from \code{glmerMod}
glmerModMultiMember <-
  setClass("glmerModMultiMember",
           contains = c("glmerMod", "merModMultiMember")
  )


#' @title Summary method for multimembership model objects
#' @param object merModMultiMember model object
#' @param ... additional arguments to be passed on to summary.merMod
#' @return summary of merModMultiMember object
#' @export
summary.merModMultiMember <- function(object, ...) {
  if (!inherits(object, "merModMultiMember")) {
    stop(
      "Cannot compute summary for objects of class: ",
      paste(class(object), collapse = ", ")
    )
  }

  # check whether model was lmer or glmer and call summary
  # if neither lmer or glmer throw error
  if (inherits(object, "lmerMod")) {
    summ <- summary(as(object, "lmerMod"), ...)
  } else {
    if (inherits(object, "glmerMod")) {
      summ <- summary(as(object, "glmerMod"), ...)
    } else {
      stop("Object class must inherit from either lmerMod or glmerMod")
    }
  }

  # store model class in summary because lme4 uses this info
  summ$objClass <- class(object)

  # add multiple membership to summary title
  summ$methTitle <- paste0(
    summ$methTitle,
    ". Model includes multiple membership random effects."
  )

  # pass multiple membership matrix through to summary
  summ$memberships <- object@memberships

  # add merModMultiMember class and return summary
  class(summ) <- c("summary.merModMultiMember", class(summ))
  return(summ)
}


#' @title Summary method for multimembership model objects
#' @param object merModMultiMember model object
#' @param ... additional arguments to be passed on to summary.merMod
#' @param ddf method for computing degrees of freedom, used by lmerTest
#' @return summary of merModMultiMember object
#' @export
summary.lmerModLmerTestMultiMember <- function(object, ..., ddf="lme4") {
  if (!inherits(object, "lmerModLmerTestMultiMember")) {
    stop(
      "Cannot compute summary for objects of class: ",
      paste(class(object), collapse = ", ")
    )
  }

  # call lmerModLmerTest summary
  if (ddf == "lme4") {
    summ <- summary(as(object, "lmerMod"), ...)
    # add merModMultiMember class
    class(summ) <- c("summary.merModMultiMember", class(summ))
  } else {
    summ <- summary(as(object, "lmerModLmerTest"), ..., ddf=ddf)
    # add lmerModMultiMemberLmerTest class
    class(summ) <- c("summary.lmerModLmerTestMultiMember", class(summ))
  }

  # store model class in summary because lme4 uses this info
  summ$objClass <- class(object)

  # add multiple membership to summary title
  summ$methTitle <- paste0(
    summ$methTitle,
    ". Model includes multiple membership random effects."
  )

  # pass multiple membership matrix through to summary
  summ$memberships <- object@memberships

  # return summary
  return(summ)
}


#' @title Print method for multiple memberships model summary
#' @param x merModMultiMember object
#' @param ... additional arguments to be passed on to print.summary.merMod
#' @export
print.summary.merModMultiMember <- function(x, ...) {
  lme4:::print.summary.merMod(x, ...)
  cat("\nGroup memberships per observation for multiple membership REs:\n")
  multimember_sums <- lapply(x$memberships, Matrix::colSums)
  print(cbind(
    "Min. per obs." = lapply(multimember_sums, min),
    "Mean per obs." = lapply(multimember_sums, mean),
    "Max. per obs." = lapply(multimember_sums, max)
  ))
  invisible(x)
}


#' @title Print method for multiple memberships model summary
#' @param x merModMultiMember object
#' @param ... additional arguments to be passed on to print.summary.merMod
#' @export
print.summary.lmerModLmerTestMultiMember <- function(x, ...) {
  print.summary.merModMultiMember(x, ...)
}


#' @title Generic for multiple dispatch of multiple membership plot functions
#' @param x object containing multimembership information, either a weight
#' matrix, multimembership model, or model summary
#' @param varname the multimembership variable to plot, if the object contains
#' more than one
#' @param ... additional arguments to pass on to Matrix::image
#' @export
plot_membership_matrix <- function(x, varname=NULL, ...) {
  UseMethod("plot_membership_matrix")
}


#' @title Plotting function for multiple membership matrix
#' @param x multimembership weight matrix
#' @param varname name of the multimembership variable, will be used as title
#' @param ... additional arguments to pass on to Matrix::image
#' @export
plot_membership_matrix.default <- function(x, varname=NULL, ...) {
  x <- Matrix::t(as(x, "dgCMatrix"))
  Matrix::image(
    x,
    main = varname,
    sub = NULL,
    ylab = "observation",
    xlab = "item",
    scales = list(
      x = list(
        at = seq(1:ncol(x)),
        labels = colnames(x),
        alternating = 3
      ),
      y = list(
        at = seq(1:nrow(x)),
        labels = rownames(x),
        alternating = 3
      )
    ),
    ...
  )
}


#' @title Plotting method for multiple memberships model
#' @param x merModMultiMember object
#' @param varname multimembership variable to be plotted
#' @param ... additional arguments to pass on to Matrix::image
#' @export
plot_membership_matrix.merModMultiMember <- function(x, varname, ...) {
  plot_membership_matrix(x@memberships[[varname]], varname, ...)
}


#' @title Plotting method for multiple memberships model summary
#' @param x merModMultiMember summary object
#' @param varname multimembership variable to be plotted
#' @param ... additional arguments to pass on to Matrix::image
#' @export
plot_membership_matrix.summary.merModMultiMember <- function(x, varname, ...) {
  plot_membership_matrix(x$memberships[[varname]], varname, ...)
}


#' @title Generic for multiple dispatch of multiple membership weights
#' histogram functions
#' @param x object containing multimembership information, either a weight
#' matrix, multimembership model, or model summary
#' @param varname the multimembership variable to plot, if the object contains
#' more than one
#' @param ... additional arguments to pass on to graphics::hist
#' @export
plot_membership_hist <- function(x, varname=NULL, ...) {
  UseMethod("plot_membership_hist")
}


#' @title Plotting function for multiple membership weights histogram
#' @param x multimembership weight matrix
#' @param varname name of the multimembership variable, will be used as title
#' @param ... additional arguments to pass on to graphics::hist
#' @export
plot_membership_hist.dgCMatrix <- function(x, varname=NULL, ...) {
  sums = Matrix::colSums(x)
  graphics::hist(
    sums,
    main = varname,
    xlab = "Number of groups associated with an observation",
    # shift breaks by .5 by default so they line up with whole integer values?
    #breaks = seq(min(sums):(max(sums) + 1)) - .5,
    ...
  )
}


#' @title Plotting method for multiple memberships model weights histogram
#' @param x merModMultiMember object
#' @param varname multimembership variable to be plotted
#' @param ... additional arguments to pass on to graphics::hist
#' @export
plot_membership_hist.merModMultiMember <- function(x, varname, ...) {
  plot_membership_hist(x@memberships[[varname]], varname, ...)
}


#' @title Plotting method for multiple memberships model summary weights
#' histogram
#' @param x merModMultiMember summary object
#' @param varname multimembership variable to be plotted
#' @param ... additional arguments to pass on to graphics::hist
#' @export
plot_membership_hist.summary.merModMultiMember <- function(x, varname, ...) {
  plot_membership_hist(x$memberships[[varname]], varname, ...)
}
