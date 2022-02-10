#' @title Multimembership random effects
#' @description Provides a wrapper around lme4::lmer to allow for multimembership random effects
#' @name lmer
#' @param formula mixed model formula
#' @param data data frame (possibly but not necessarily containing factors
#' @param memberships list of weights matrices  with which to replace Zt components
#' @param ... additional arguments to pass through to lme4::lmer
#' @return lme4 model object
#' @export
#' @import lme4
#' @examples
#' df <- data.frame(
#'   x = seq(60) + runif(60, 0, 10),
#'   y = seq(60) + rep(runif(6, 0, 10), 10),
#'   memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
#' )
#' weights <- weights_from_vector(df$memberships)
#' lmer(y ~ x + (1|members), df, memberships=list(members=weights))
lmer <- function(formula, data = NULL, REML = TRUE, control = lme4::lmerControl(),
                 start = NULL, verbose = 0L, weights = NULL, na.action=na.omit,
                 offset = NULL, contrasts = NULL, devFunOnly = FALSE, memberships=NULL) {
  ## TODO: call weights_from_vector from within this function
  ## FIXME: docstring

  orig_call <- match.call()
  # TODO: detect if lmerTest is loaded
  if(is.null(memberships)){
    return(lme4::lmer(formula, data=data, REML=REML, control=control,
                      start=start, verbose=verbose, weights=weights,
                      na.action=na.action, offset=offset, contrasts=contrasts,
                      devFunOnly=devFunOnly))
  }

  ## FIXME: pass ... through appropriately -> not all args should be passed to same function?
  # so: use lme4::lmer to figure out which args need to be passed through, and where they go
  # e.g. control and start go to optimizerLmer()
  # but devFunOnly needs to be caught in an if/else block
  # and where does the arg contrasts go? etc. etc.

  ## FIXME: test dimensions

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
    w <- which(mnms==gvars[i])

    if (length(w) > 0) {

      M <- Matrix::Matrix(memberships[[w]])

      ## extract LHS (effect)
      form <- as.formula(substitute(~z, list(z=fb[[i]][[2]])))

      ## construct model matrix & compute Khatri-Rao product
      X <- model.matrix(form, data=data)
      Zt <- Matrix::KhatriRao(M, t(X), make.dimnames=TRUE)

      ## FIXME: mess with names?
      Ztlist[[fbnm]] <- Zt

      ## if necessary, add factor to data
      if (!gvars[i] %in% names(data)) {
        ## if the factor has non-trivial ordering, it should be included
        ## in the data.  Do we have to worry about ordering of Z? test!
        data[[gvars[i]]] <- rep_len(factor(rownames(memberships[[w]])), dim(data)[1])
      }
    }
  }

  lmod <- lFormula(formula, data=data, REML=REML,
                   weights=weights, na.action=na.action, offset=offset,
                   contrasts=contrasts, control=control)

  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)

  ## finish fitting
  devfun <- do.call(mkLmerDevfun, lmod)
  if(devFunOnly) return(devfun)
  opt <- optimizeLmer(devfun,
                      optimizer    = control$optimizer,
                      restart_edge = control$restart_edge,
                      boundary.tol = control$boundary.tol,
                      start = NULL,
                      verbose = 0L,
                      control = control$optCtrl)

  m1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr=lmod$fr, mc=orig_call)

  # convert model object to lmerModMultiMember -> does this need any additional slots/attributes?
  # or do we just want to output some additional information when summary.lmerModMultiMember is called?
  # e.g. for each random effect: don't show number of levels, but number of groups & number of unique group members
  m1 <- as(m1, "lmerModMultiMember")
  return(m1)
}

#' @importFrom lme4 fixef
#' @export
lme4::fixef

#' @importFrom lme4 ranef
#' @export
lme4::ranef


#' @title Model object for multimembership mixed models
#' @description The \code{lmerModMultiMember} class extends \code{lmerMod} (which extends
#' \code{merMod}) from the \pkg{lme4}-package.
#' @seealso \code{\link[lme4]{lmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @importClassesFrom lme4 lmerMod
#' @return An object of class \code{lmerModMultiMember} similar to
#' \code{lmerMod} objects (see \code{\link[lme4]{merMod}}) with extra information
#' about the multimembership random effects.
lmerModMultiMember <-
  setClass("lmerModMultiMember",
           contains = c("lmerMod")
           )
