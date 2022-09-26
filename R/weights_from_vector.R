#' @title Weight matrix from vector of memberships
#' @description Provides a helper function for generating weight matrices
#' @name weights_from_vector
#' @param membership_vector vector of strings containing group memberships
#' @param sep separator that delimits group memberships in membership_vector
#' @return sparse weight/indicator matrix of type Matrix::dgCmatrix
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
#' @return sparse weight/indicator matrix of type Matrix::dgCmatrix
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
#' @param a sparse weight/indicator matrix of class Matrix::dgCmatrix
#' @param b sparse weight/indicator matrix of class Matrix::dgCmatrix
#' @return sparse interaction weight/indicator matrix of class Matrix::dgCmatrix
#' @export
#' @examples
#' a <- rep(c("k", "l", "l,k"), 2)
#' b <- rep(c("m", "n"), 3)
#' Wa <- weights_from_vector(a)
#' Wb <- Matrix::fac2sparse(b)
#' Wab <- interaction_weights(Wa, Wb)
interaction_weights <- function(a, b) {
  abrows <- as.character(interaction(expand.grid(rownames(a), rownames(b))))
  ab <- Matrix::Matrix(0, nrow = length(abrows), ncol = ncol(a),
                       dimnames = list(abrows, colnames(a)))
  for (i in seq_along(rownames(a))) {
    for (j in seq_along(rownames(b))) {
      ab[((j - 1) * length(rownames(a))) + i, ] <-
        a[i, , drop = FALSE] * b[j, , drop = FALSE]
    }
  }
  return(ab)
}


idx_to_weights <- function(idx, nobs) {

  # get sorted vector of unique groups
  groups <- sort(unique(idx$values))

  # convert group membership indices to numeric factors
  idx[] <- lapply(idx, function(x) as.numeric(factor(x)))

  # initialize sparse weight matrix with zeros
  weights <- Matrix::Matrix(
    0,
    nrow = length(groups),
    ncol = nobs,
    dimnames = list(as.character(groups), as.character(seq(nobs)))
  )

  # fill in sparse weight matrix with 1s for each group membership
  weights[as.matrix(idx)] <- 1

  # return sparse weight matrix
  return(weights)
}
