#' @title Weight matrix from vector of memberships
#' @description Provides a helper function for generating weight matrices
#' @name weights_from_vector
#' @param membership_vector vector of strings containing group memberships
#' @param sep separator that delimits group memberships in membership_vector
#' @return sparse weight/indicator matrix of type Matrix::dgCMatrix
#' @export
#' @import stats utils methods
#' @examples
#' a <- c("k,l,m", "k,m", "k", "l", "l,k")
#' Wa <- weights_from_vector(a)
weights_from_vector <- function(membership_vector, sep = ",") {

  # get number of observations
  nobs <- length(membership_vector)

  # split membership vector into list of vectors
  membership_list <- strsplit(membership_vector, split = sep)

  # create list of group membership indices
  idx <- stack(setNames(membership_list, seq_along(membership_list)))

  # hand off to idx_to_weights
  return(idx_to_weights(idx, nobs))
}


#' @title Weight matrix from columns containing memberships
#' @description Provides a helper function for generating weight matrices
#' @name weights_from_columns
#' @param membership_columns columns containing group memberships
#' @return sparse weight/indicator matrix of type Matrix::dgCMatrix
#' @export
#' @examples
#' a <- cbind(
#'   c("k", "k", "l", NA, "m"),
#'   c("l", "m", "k", "m", "k"),
#'   c("m", NA, "m", NA, NA)
#' )
#' Wa <- weights_from_columns(a)
weights_from_columns <- function(membership_columns) {

  # get number of observations
  nobs <- nrow(membership_columns)

  # create list of group membership indices
  idx <- na.omit(stack(setNames(
    as.matrix(membership_columns),
    rep(1:nobs, ncol(membership_columns))
  )))

  # hand off to idx_to_weights
  return(idx_to_weights(idx, nobs))
}


#' @title Weight matrix from the interaction of two other weight matrices
#' @description Provides a helper function for pre-generating interactions.
#' Takes sparse weight/indicator matrices generated with e.g.
#' Matrix::fac2sparse() or weights_from_vector() as input.
#' @name interaction_weights
#' @param a sparse weight/indicator matrix of class Matrix::dgCMatrix
#' @param b sparse weight/indicator matrix of class Matrix::dgCMatrix
#' @return sparse interaction weight/indicator matrix of class Matrix::dgCMatrix
#' @export
#' @examples
#' a <- rep(c("k", "l", "l,k"), 2)
#' b <- rep(c("m", "n"), 3)
#' Wa <- weights_from_vector(a)
#' Wb <- Matrix::fac2sparse(b)
#' Wab <- interaction_weights(Wa, Wb)
interaction_weights <- function(a, b) {

  # check if a and b matrices contain same number of cases
  if (ncol(a) != ncol(b)) {
    stop("Matrices must have same number of observations/cases.")
  }

  # expand a and b rownames into rows for an ab matrix
  abrows <- as.character(interaction(expand.grid(rownames(a), rownames(b))))

  # expand a and b matrices by repeating their rows to match ab rows
  # note a repeats using times pattern, while b repeats using each pattern
  # then multiply a and b together to form ab
  ab <- a[rep(1:nrow(a), times=length(abrows) / nrow(a)), ] *
    b[rep(1:nrow(b), each=length(abrows) / nrow(b)), ]

  # assign ab rownames
  rownames(ab) <- abrows

  # return ab matrix, with rows sorted alphabetically by row name
  # this matches lme4 internal behavior
  return(ab[sort(abrows), ])
}


#' @title Indicator matrix from two vectors for Bradley-Terry models
#' @description Provides a helper function for generating indicator matrices
#' @name bradleyterry_from_vectors
#' @param winners vector of strings containing the winner for each observation
#' @param losers vector of strings containing the loser for each observation
#' @return sparse Bradley-Terry indicator matrix of type Matrix::dgCMatrix
#' @export
#' @examples
#' winners <- c("k", "k", "l", "m", "m")
#' losers <- c("l", "m", "m", "k", "l")
#' Wwl <- bradleyterry_from_vectors(winners, losers)
bradleyterry_from_vectors <- function(winners, losers) {
  return(bradleyterry_from_sparse(Matrix::fac2sparse(winners),
                                  Matrix::fac2sparse(losers)))
}


#' @title Indicator matrix from two sparse matrices for Bradley-Terry models
#' @description Provides a helper function for generating indicator matrices
#' @name bradleyterry_from_sparse
#' @param winners sparse matrix containing the winner for each observation
#' @param losers sparse matrix containing the loser for each observation
#' @return sparse Bradley-Terry indicator matrix of type Matrix::dgCMatrix
#' @export
#' @examples
#' winners <- Matrix::fac2sparse(c("k", "k", "l", "m", "m"))
#' losers <- Matrix::fac2sparse(c("l", "m", "m", "k", "l"))
#' Wwl <- bradleyterry_from_sparse(winners, losers)
bradleyterry_from_sparse <- function(winners, losers) {

  # reassign to shorter variable names
  Wa <- winners
  Wb <- losers

  # get number of observations
  nobs <- ncol((Wa))

  # create empty sparse matrix with rows from union of Wa and Wb rownames
  Wab_rownames <- sort(union(rownames(Wa), rownames(Wb)))
  Wab <- Matrix::Matrix(
    0,
    nrow = length(Wab_rownames),
    ncol = ncol(Wa),
    dimnames = list(Wab_rownames, as.character(1:nobs))
  )

  # fill sparse matrix with Wa indicators and negative Wb indicators
  Wab[rownames(Wa), ] <- Wa
  Wab[rownames(Wb), ] <- Wab[rownames(Wb), ] - Wb

  # return sparse Bradley-Terry indicator matrix
  return(Wab)
}


idx_to_weights <- function(idx, nobs) {

  # get sorted vector of unique groups
  groups <- sort(unique(idx$values))

  # convert group membership indices to numeric factors
  idx[] <- lapply(idx, function(x) as.numeric(factor(x)))

  # sum values for each index
  # (e.g. for when an observation has multiple members of the same group/class)
  idx$count <- 1
  idx <- aggregate(count ~ values + ind, data = idx, FUN = sum)

  # initialize sparse weight matrix with zeros
  weights <- Matrix::Matrix(
    0,
    nrow = length(groups),
    ncol = nobs,
    dimnames = list(as.character(groups), as.character(1:nobs))
  )

  # fill in sparse weight matrix with 1s for each group membership
  weights[as.matrix(idx[, c("values", "ind")])] <- idx$count

  # return sparse weight matrix
  return(weights)
}
