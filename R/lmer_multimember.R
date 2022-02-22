#' @title multimembership linear mixed effects models
#' @description lme4::lmer but with multimembership random effects
#' @name lmer
#' @inheritParams lme4::lmer
#' @param memberships list of weights matrices to replace Zt components with
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

  # TODO: detect if lmerTest is loaded

  if (is.null(memberships)) {
    return(lme4::lmer(formula,
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
    ))
  }

  # get names of multimembership variables
  mnms <- names(memberships)

  # get index of bars (location of random effects) in model formula
  fb <- lme4::findbars(formula)

  # get random effect grouping variables
  gvars <- vapply(fb, function(x) deparse(x[[3]]), character(1))

  Ztlist <- list()

  # iterate over random effects
  for (i in seq_along(fb)) {
    fbnm <- deparse(fb[[i]])

    ## find corresponding random-effects term
    w <- which(mnms == gvars[i])

    if (length(w) > 0) {
      M <- Matrix::Matrix(memberships[[w]])

      ## extract LHS (effect)
      form <- as.formula(substitute(~z, list(z = fb[[i]][[2]])))

      ## construct model matrix & compute Khatri-Rao product
      X <- model.matrix(form, data = data)
      Zt <- Matrix::KhatriRao(M, t(X), make.dimnames = TRUE)

      ## FIXME: mess with names?
      Ztlist[[fbnm]] <- Zt

      ## if necessary, add factor to data
      if (!gvars[i] %in% names(data)) {
        ## if the factor has non-trivial ordering, it should be included
        ## in the data.  Do we have to worry about ordering of Z? test!
        data[[gvars[i]]] <- rep_len(
          factor(rownames(memberships[[w]])),
          dim(data)[1]
        )
      }
    }
  }

  lmod <- lFormula(formula,
    data = data,
    REML = REML,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    control = control
  )

  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)

  ## finish fitting
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

  m1 <- mkMerMod(environment(devfun),
    opt,
    lmod$reTrms,
    fr = lmod$fr,
    mc = orig_call
  )

  # convert model object to lmerModMultiMember
  # TODO: does this object need any additional slots/attributes?
  # or do we just want to output some additional information when
  # summary.lmerModMultiMember is called?
  # e.g. for each random effect: don't show number of levels,
  # but number of groups & number of unique group members
  m1 <- as(m1, "lmerModMultiMember")
  m1@memberships <- memberships
  return(m1)
}


#' @title multimembership generalized linear mixed effects models
#' @description lme4::glmer but with multimembership random effects
#' @name glmer
#' @inheritParams lme4::glmer
#' @param memberships list of weights matrices to replace Zt components with
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

  # TODO: detect if lmerTest is loaded

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
  mnms <- names(memberships)

  # get index of bars (location of random effects) in model formula
  fb <- lme4::findbars(formula)

  # get random effect grouping variables
  gvars <- vapply(fb, function(x) deparse(x[[3]]), character(1))

  Ztlist <- list()

  # iterate over random effects
  for (i in seq_along(fb)) {
    fbnm <- deparse(fb[[i]])

    ## find corresponding random-effects term
    w <- which(mnms == gvars[i])

    if (length(w) > 0) {
      M <- Matrix::Matrix(memberships[[w]])

      ## extract LHS (effect)
      form <- as.formula(substitute(~z, list(z = fb[[i]][[2]])))

      ## construct model matrix & compute Khatri-Rao product
      X <- model.matrix(form, data = data)
      Zt <- Matrix::KhatriRao(M, t(X), make.dimnames = TRUE)

      ## FIXME: mess with names?
      Ztlist[[fbnm]] <- Zt

      ## if necessary, add factor to data
      if (!gvars[i] %in% names(data)) {
        ## if the factor has non-trivial ordering, it should be included
        ## in the data.  Do we have to worry about ordering of Z? test!
        data[[gvars[i]]] <- rep_len(
          factor(rownames(memberships[[w]])),
          dim(data)[1]
        )
      }
    }
  }

  glmod <- glFormula(formula,
    data = data,
    family = family,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    control = control
  )

  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    glmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  glmod$reTrms$Zt <- do.call(rbind, glmod$reTrms$Ztlist)

  ## finish fitting
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

  ## FIX ME: allow calc.derivs, use.last.params etc. if nAGQ=0
  if (control$nAGQ0initStep) {
    opt <- optimizeGlmer(devfun,
      optimizer = control$optimizer[[1]],
      ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
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


    ## update deviance function to include fixed effects as inputs
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
    ## reoptimize deviance function over covariance parameters and fixed effects
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

  # convert model object to glmerModMultiMember
  # does this need any additional slots/attributes?
  # or do we just want to output some additional information when
  # summary.glmerModMultiMember is called?
  # e.g. for each random effect: don't show number of levels,
  # but number of groups & number of unique group members
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
#' @description The \code{lmerModMultiMember} class extends \code{lmerMod}
#' (which extends \code{merMod}) from the \pkg{lme4}-package.
#' @seealso \code{\link[lme4]{lmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @importClassesFrom lme4 lmerMod
#' @return An object of class \code{lmerModMultiMember} similar to
#' \code{lmerMod} objects (see \code{\link[lme4]{merMod}}) with
#' extra information about the multimembership random effects.
lmerModMultiMember <-
  setClass("lmerModMultiMember",
    contains = c("lmerMod"), slots = c(memberships = "list")
  )


#' @title Model object for multimembership generalized linear mixed models
#' @description The \code{glmerModMultiMember} class extends \code{glmerMod}
#' (which extends \code{merMod}) from the \pkg{lme4}-package.
#' @seealso \code{\link[lme4]{glmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @importClassesFrom lme4 glmerMod
#' @return An object of class \code{glmerModMultiMember} similar to
#' \code{glmerMod} objects (see \code{\link[lme4]{merMod}}) with
#' extra information about the multimembership random effects.
glmerModMultiMember <-
  setClass("glmerModMultiMember",
    contains = c("glmerMod"), slots = c(memberships = "list")
  )


#' @title summary method for multimembership model objects
#' @export
summary.lmerModMultiMember <- function(object, ...) {
  # ddf <- match.arg(ddf)
  if(!inherits(object, "lmerModMultiMember") && !inherits(object, "lmerMod")) {
    stop("Cannot compute summary for objects of class: ",
         paste(class(object), collapse = ", "))
  }

  ##############################################################################
  ## lmerTest block, unclear if we really need this:
  ## what is the point of coercing to lmerModMultiMember if the model wasn't
  ## a multimembership model in the first place?
  ##############################################################################
  #if(!inherits(object, "lmerModMultiMember") && inherits(object, "lmerMod")) {
  #  message("Coercing object to class 'lmerModMultiMember'")
  #  object <- as_lmerModMultiMember(object)
  #  if(!inherits(object, "lmerModMultiMember")) {
  #    warning("Failed to coerce object to class 'lmerModMultiMember'")
  #    return(summary(object))
  #  }
  #}
  ##############################################################################
  memberships <- object@memberships
  summ <- summary(as(object, "lmerMod"), ...)
  summ$objClass <- class(object) # Used by lme4:::print.summary.lmerMod
  summ$methTitle <- paste0(summ$methTitle,
                           " with multiple membership random effects")
  summ$memberships <- memberships
  class(summ) <- c("summary.lmerModMultiMember", class(summ))

  summ
}
