#' @title Multimembership random effects
#' @description Provides a wrapper around lme4::lmer to allow for multimembership random effects
#' @name lmer_multimember
#' @param formula mixed model formula
#' @param data data frame (possibly but not necessarily containing factors
#' @param memberships list of weights matrices  with which to replace Zt components
#' @return lme4 model object
#' @export
#' @examples
#' df <- data.frame(
#'   x = seq(60) + runif(60, 0, 10),
#'   y = seq(60) + rep(runif(6, 0, 10), 10),
#'   memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
#' )
#' weights <- weights_from_vector(df$memberships)
#' lmer_multimember(y ~ x + (1|members), df, memberships=list(members=weights))
lmer_multimember <- function(formula, data=NULL, REML = TRUE,
                 control = lme4::lmerControl(), start = NULL
                 , verbose = 0L
                 , subset, weights, na.action, offset
                 , contrasts = NULL
                 , devFunOnly=FALSE
                 , memberships = list()
) {
  mc <- mcout <- match.call()
  missCtrl <- missing(control)
  ## see functions in modular.R for the body ..
  if (!missCtrl && !inherits(control, "lmerControl")) {
    if(!is.list(control)) stop("'control' is not a list; use lmerControl()")
    ## back-compatibility kluge
    warning("passing control as list is deprecated: please use lmerControl() instead",
            immediate.=TRUE)
    control <- do.call(lme4::lmerControl, control)
  }
  mc$control <- control ## update for  back-compatibility kluge

  ################ START lmerMultiMember block ################
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
  ################ END lmerMultiMember block ################

  ## https://github.com/lme4/lme4/issues/50
  ## parse data and formula
  #mc[[1]] <- quote(lme4::lFormula)
  #lmod <- eval(mc, parent.frame(1L))
  #mcout$formula <- lmod$formula
  #lmod$formula <- NULL
  # not sure why this parent.frame construction was used here
  # but it messes up the added (fake) factors that we need for the multimembership model to work right

  ################ START lmerMultiMember block ################
  ## parse the formula and data
  lmod <- lme4::lFormula(formula, data)
  mcout$formula <- lmod$formula
  lmod$formula <- NULL

  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)
  ################ END lmerMultiMember block ################

  ## create deviance function for covariance parameters (theta)
  devfun <- do.call(lme4::mkLmerDevfun,
                    c(lmod,
                      list(start=start, verbose=verbose, control=control)))
  if (devFunOnly) return(devfun)
  ## optimize deviance function over covariance parameters
  if (identical(control$optimizer,"none"))
    stop("deprecated use of optimizer=='none'; use NULL instead")
  opt <- if (length(control$optimizer)==0) {
    s <- getStart(start, environment(devfun)$pp)
    list(par=s,fval=devfun(s),
         conv=1000,message="no optimization")
  }  else {
    lme4::optimizeLmer(devfun, optimizer = control$optimizer,
                 restart_edge = control$restart_edge,
                 boundary.tol = control$boundary.tol,
                 control = control$optCtrl,
                 verbose=verbose,
                 start=start,
                 calc.derivs=control$calc.derivs,
                 use.last.params=control$use.last.params)
  }
  cc <- lme4::checkConv(attr(opt,"derivs"), opt$par,
                  ctrl = control$checkConv,
                  lbound = environment(devfun)$lower)

  ## prepare output
  return(as(lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr,
           mc = mcout, lme4conv=cc), "lmerModMultiMember"))
}

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
