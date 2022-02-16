#' @title Weight matrix from vector of memberships
#' @description Provides a helper function for generating weight matrices
#' @name weights_from_vector
#' @param membership_vector vector of strings containing group memberships
#' @param sep separator that delimits group memberships in membership_vector
#' @return sparse weight matrix
#' @export
#' @examples
#' member_vec <-  c("a,b,c", "a,c", "a", "b", "b,a")
#' weights_from_vector(member_vec)
weights_from_vector <- function(membership_vector, sep=",") {

  # get number of observations
  nobs <- length(membership_vector)

  # split membership vector into list of vectors
  membership_list <- strsplit(membership_vector, split=sep)

  # create list of group membership indices
  idx <- stack(setNames(membership_list, seq_along(membership_list)))

  # hand off to idx_to_weights
  return(idx_to_weights(idx, nobs))
}


#' @title Weight matrix from columns containing memberships
#' @description Provides a helper function for generating weight matrices
#' @name weights_from_columns
#' @param membership_columns columns containing group memberships
#' @return sparse weight matrix
#' @export
#' @examples
#' member_cols <- cbind(
#'   c("a", "a", "b", NA, "c"),
#'   c("b", "c", "a", "c", "a"),
#'   c("c", NA, "c", NA, NA)
#' )
#' weights_from_columns(member_cols)
weights_from_columns <- function(membership_columns) {

  # get number of observations
  nobs <- nrow(membership_columns)

  # create list of group membership indices
  idx <- na.omit(stack(setNames(membership_columns,
                                rep(1:nobs, ncol(membership_columns)))))

  # hand off to idx_to_weights
  return(idx_to_weights(idx, nobs))
}


idx_to_weights <- function(idx, nobs) {

  # get sorted vector of unique groups
  groups <- sort(unique(idx$values))

  # convert group membership indices to numeric factors
  idx[] <- lapply(idx, function(x) as.numeric(factor(x)))

  # initialize sparse weight matrix with zeros
  weights <- Matrix::Matrix(
    0,
    nrow=length(groups),
    ncol=nobs,
    dimnames=list(as.character(groups), as.character(seq(nobs)))
  )

  # fill in sparse weight matrix with 1s for each group membership
  weights[as.matrix(idx)] <- 1

  # return sparse weight matrix
  return(weights)
}
