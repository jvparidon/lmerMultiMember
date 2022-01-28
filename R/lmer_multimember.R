#' @title Multimembership random effects
#' @description Provides a wrapper around lme4::lmer to allow for multimembership random effects
#' @name lmer_multimember
#' @aliases lmer_multimember
#' @param formula mixed model formula
#' @param data data frame (possibly but not necessarily containing factors
#' @param memb_mat list of weights matrices  with which to replace Zt components
#' @param ... additional arguments to pass through to lme4::lmer
#' @return lme4 model object
#' @examples
#' df <- data.frame(
#'   x = seq(5),
#'   y = seq(5),
#'   memberships = c("a,b,c", "a,c", "a", "b", "b,a")
#' )
#' weights <- weights_from_vector(df$memberships)
#' lmer_multimember(y ~ x + (1|members), df, memb_mat=list(members=weights))
lmer_multimember <- function(formula, data, memb_mat=list(), ...) {
  ## FIXME: pass ... through appropriately
  ## FIXME: test dimensions
  # get names of multimembership variables
  mnms <- names(memb_mat)
  # get index of bars (location of random effects) in model formula
  fb <- findbars(formula)
  # get random effect grouping variables
  gvars <- vapply(fb, function(x) deparse(x[[3]]), character(1))
  Ztlist <- list()
  for (i in seq_along(fb)) {
    fbnm <- deparse(fb[[i]])
    ## find corresponding random-effects term
    w <- which(mnms==gvars[i])
    if (length(w)>0) {
      M <- Matrix::Matrix(memb_mat[[w]])
      ## extract LHS (effect)
      form <- as.formula(substitute(~z,list(z=fb[[i]][[2]])))
      ## construct model matrix & compute Khatri-Rao product
      X <- model.matrix(form,data=data)
      Zt <- Matrix::KhatriRao(t(M),t(X),make.dimnames=TRUE)
      ## FIXME: mess with names?
      Ztlist[[fbnm]] <- Zt
      ## if necessary, add factor to data
      if (!gvars[i] %in% names(data)) {
        ## if the factor has non-trivial ordering, it should be included
        ## in the data.  Do we have to worry about ordering of Z? test!
        data[[gvars[i]]] <- rep_len(factor(colnames(memb_mat[[w]])), dim(data)[1])
      }
    } ## if  (length(w)>0)
  } ## for i in seq(fb)
  ## call lFormula  (FIXME: allow glFormula)
  lmod <- lFormula(formula,data=data)
  ## substitute new Ztlist elements
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind,lmod$reTrms$Ztlist)
  ## finish fitting
  devfun <- do.call(mkLmerDevfun, lmod)
  opt <- optimizeLmer(devfun)
  m1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
  return(m1)
}

#' @title Weight matrix from vector of memberships
#' @description Provides a helper function for generating weight matrices
#' @name weights_from_vector
#' @aliases weights_from_vector
#' @param membership_vector vector of strings containing group memberships
#' @param sep separator that delimits group memberships in membership_vector
#' @return sparse weight matrix
#' @examples
#' member_vec <-  c("a,b,c", "a,c", "a", "b", "b,a")
#' weights_from_vector(member_vec)
weights_from_vector <- function(membership_vector, sep=",") {
  # split membership vector into list of vectors
  # (is this necessary? could also just require a list of vectors as input to begin with)
  membership_list <- strsplit(membership_vector, split=sep)
  # create list of group membership indices
  idx <- stack(setNames(membership_list, seq_along(membership_list)))
  # get sorted vector of unique groups
  groups <- sort(unique(idx$values))
  # convert group membership indices to numeric factors
  idx[] <- lapply(idx, function(x) as.numeric(factor(x)))
  # initialize sparse weight matrix with zeros
  weights <- Matrix(
    0,
    nrow=length(groups),
    ncol=length(membership_list),
    dimnames=list(as.character(groups), as.character(seq(length(membership_list))))
  )
  # fill in sparse weight matrix with 1s for each group membership
  weights[as.matrix(idx)] <- 1
  # return transposed sparse weight matrix
  return(t(weights))
}
