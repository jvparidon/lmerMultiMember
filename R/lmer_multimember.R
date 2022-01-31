#' @title Linear mixed effects model with multimembership random effects
#' @description Version of lme4::lmer to allow for multimembership random effects
#' @name lmer_multimember
#' @param formula mixed model formula
#' @param data data frame (possibly but not necessarily containing factors
#' @param memb_mat list of weights matrices  with which to replace Zt components
#' @param ... additional arguments to pass through to lme4::lmer
#' @return lme4 model object
#' @export
#' @examples
#' df <- data.frame(
#'   x = seq(5),
#'   y = seq(5),
#'   memberships = c("a,b,c", "a,c", "a", "b", "b,a")
#' )
#' weights <- weights_from_vector(df$memberships)
#' lmer_multimember(y ~ x + (1|members), df, memb_mat=list(members=weights))
lmer_multimember <- function(formula, data, memb_mat=list(), ...) {

  ## FIXME: pass ... through appropriately -> not all args should be passed to same function?
  # so: use lme4::lmer to figure out which args need to be passed through, and where they go
  # e.g. control and start go to optimizerLmer()
  # but devFunOnly needs to be caught in an if/else block
  # and where does the arg contrasts go? etc. etc.

  ## FIXME: test dimensions

  # get names of multimembership variables
  mnms <- names(memb_mat)

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

      M <- Matrix::Matrix(memb_mat[[w]])

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
        data[[gvars[i]]] <- rep_len(factor(rownames(memb_mat[[w]])), dim(data)[1])
      }
    }
  }

  ## call lFormula  (FIXME: allow glFormula) -> is it easiest to have a separate lmerMultiMember::glmer function
  # that passes all its arguments through to lmerMultiMember::lmer, and then catch that here
  # and run it through lme4::glFormula instead? (will that work?)
  lmod <- lme4::lFormula(formula, data=data)

  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind, lmod$reTrms$Ztlist)

  ## finish fitting
  devfun <- do.call(lme4::mkLmerDevfun, lmod)
  opt <- lme4::optimizeLmer(devfun)
  m1 <- lme4::mkMerMod(environment(devfun), opt, lmod$reTrms, fr=lmod$fr)

  # convert model object to lmerModMultiMember -> does this need any additional slots/attributes?
  # or do we just want to output some additional information when summary.lmerModMultiMember is called?
  # e.g. for each random effect: don't show number of levels, but number of groups & number of unique group members
  m1 <- as(m1, "lmerModMultiMember")
  return(m1)
}


#' @title Generalized linear mixed effects model with multimembership random effects
#' @description Version of lme4::glmer that allow for multimembership random effects
#' @name glmer_multimember
#' @param formula mixed model formula
#' @param data data frame (possibly but not necessarily containing factors
#' @param memberships list of weights matrices  with which to replace Zt components
#' @return lme4 model object
#' @export
#' @importFrom lme4 GHrule
#' @examples
#' df <- data.frame(
#'   x = runif(60, 0, 1),
#'   y = rbinom(60, 1, 0.6),
#'   memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
#' )
#' weights <- weights_from_vector(df$memberships)
#' glmer_multimember(y ~ x + (1|members), df, family=binomial, memberships=list(members=weights))
glmer_multimember <- function(formula, data=NULL
                  , family = gaussian
                  , control = lme4::glmerControl()
                  , start = NULL
                  , verbose = 0L
                  , nAGQ = 1L
                  , subset, weights, na.action, offset, contrasts = NULL
                  , mustart, etastart
                  , devFunOnly = FALSE
                  , memberships = list())
{
  if (!inherits(control, "glmerControl")) {
    if(!is.list(control)) stop("'control' is not a list; use glmerControl()")
    ## back-compatibility kluge
    if (class(control)[1]=="lmerControl") {
      warning("please use glmerControl() instead of lmerControl()",
              immediate.=TRUE)
      control <-
        ## unpack sub-lists
        c(control[!names(control) %in% c("checkConv","checkControl")],
          control$checkControl,control$checkConv)
      control["restart_edge"] <- NULL ## not implemented for glmer
    } else {
      msg <- "Use control=glmerControl(..) instead of passing a list"
      if(length(cl <- class(control))) {
        msg <- paste(msg, "of class", dQuote(cl[1]))
      }
      warning(msg, immediate.=TRUE)
    }

    control <- do.call(lme4::glmerControl, control)
  }
  mc <- mcout <- match.call()

  ## family-checking code duplicated here and in glFormula (for now) since
  ## we really need to redirect at this point; eventually deprecate formally
  ## and clean up
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame(2))
  if( is.function(family)) family <- family()
  if (isTRUE(all.equal(family, gaussian()))) {
    ## redirect to lmer (with warning)
    warning("calling glmer() with family=gaussian (identity link) as a shortcut to lmer() is deprecated;",
            " please call lmer() directly")
    mc[[1]] <- quote(lmer_multimember)
    mc["family"] <- NULL            # to avoid an infinite loop
    return(eval(mc, parent.frame()))
  }

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

  ## parse the formula and data
  #mc[[1]] <- quote(lme4::glFormula)
  #glmod <- eval(mc, parent.frame(1L))
  #mcout$formula <- glmod$formula
  #glmod$formula <- NULL
  # not sure why this parent.frame construction was used here
  # but it messes up the added (fake) factors that we need for the multimembership model to work right

  ################ START lmerMultiMember block ################
  ## parse the formula and data
  glmod <- lme4::glFormula(formula, data, family = binomial)
  mcout$formula <- glmod$formula
  glmod$formula <- NULL

  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    glmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  glmod$reTrms$Zt <- do.call(rbind, glmod$reTrms$Ztlist)
  ################ END lmerMultiMember block ################

  ## create deviance function for covariance parameters (theta)

  nAGQinit <- if(control$nAGQ0initStep) 0L else 1L
  devfun <- do.call(lme4::mkGlmerDevfun, c(glmod, list(verbose = verbose,
                                                 control = control,
                                                 nAGQ = nAGQinit)))
  if (nAGQ==0 && devFunOnly) return(devfun)
  ## optimize deviance function over covariance parameters

  if (is.list(start)) {
    start.bad <- setdiff(names(start),c("theta","fixef"))
    if (length(start.bad)>0) {
      stop(sprintf("bad name(s) for start vector (%s); should be %s and/or %s",
                   paste(start.bad,collapse=", "),
                   shQuote("theta"),
                   shQuote("fixef")),call.=FALSE)
    }
    if (!is.null(start$fixef) && nAGQ==0)
      stop("should not specify both start$fixef and nAGQ==0")
  }

  ## FIX ME: allow calc.derivs, use.last.params etc. if nAGQ=0
  if(control$nAGQ0initStep) {
    opt <- lme4::optimizeGlmer(devfun,
                         optimizer = control$optimizer[[1]],
                         ## DON'T try fancy edge tricks unless nAGQ=0 explicitly set
                         restart_edge=if (nAGQ==0) control$restart_edge else FALSE,
                         boundary.tol=if (nAGQ==0) control$boundary.tol else 0,
                         control = control$optCtrl,
                         start=start,
                         nAGQ = 0,
                         verbose=verbose,
                         calc.derivs=FALSE)
  }

  if(nAGQ > 0L) {


    ## update deviance function to include fixed effects as inputs
    devfun <- lme4::updateGlmerDevfun(devfun, glmod$reTrms, nAGQ = nAGQ)

    ################ END lmerMultiMember block ################
    # updateStart is a utility function from lme4 that is necessary here
    # but for some reason lme4 doesn't export it
    updateStart <- function(start, theta) {
      if (is.numeric(start)) {
        theta
      } else {
        if (!is.null(start$theta))
          start$theta <- theta
        start
      }
    }
    ################ END lmerMultiMember block ################

    if (control$nAGQ0initStep) {
      start <- updateStart(start,theta=opt$par)
    }
    ## if nAGQ0 was skipped
    ## we don't actually need to do anything here, it seems --
    ## getStart gets called again in optimizeGlmer

    if (devFunOnly) return(devfun)
    ## reoptimize deviance function over covariance parameters and fixed effects
    opt <- lme4::optimizeGlmer(devfun,
                         optimizer = control$optimizer[[2]],
                         restart_edge=control$restart_edge,
                         boundary.tol=control$boundary.tol,
                         control = control$optCtrl,
                         start=start,
                         nAGQ=nAGQ,
                         verbose = verbose,
                         stage=2,
                         calc.derivs=control$calc.derivs,
                         use.last.params=control$use.last.params)
  }
  cc <- if (!control$calc.derivs) NULL else {
    if (verbose > 10) cat("checking convergence\n")
    lme4::checkConv(attr(opt,"derivs"),opt$par,
              ctrl = control$checkConv,
              lbound=environment(devfun)$lower)
  }

  ## prepare output
  return(as(lme4::mkMerMod(environment(devfun), opt, glmod$reTrms, fr = glmod$fr,
           mc = mcout, lme4conv=cc), "glmerModMultiMember"))

}

#' @title Model object for multimembership linear mixed models
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

#' @title Model object for multimembership generalized linear mixed models
#' @description The \code{glmerModMultiMember} class extends \code{glmerMod} (which extends
#' \code{merMod}) from the \pkg{lme4}-package.
#' @seealso \code{\link[lme4]{glmer}} and \code{\link[lme4]{merMod}}
#' @export
#' @importClassesFrom lme4 glmerMod
#' @return An object of class \code{glmerModMultiMember} similar to
#' \code{glmerMod} objects (see \code{\link[lme4]{merMod}}) with extra information
#' about the multimembership random effects.
glmerModMultiMember <-
  setClass("glmerModMultiMember",
           contains = c("glmerMod")
  )
