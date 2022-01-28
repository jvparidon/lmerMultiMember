#' @title Weight matrix from vector of memberships
#' @description Provides a helper function for generating weight matrices
#' @name weights_from_vector
#' @aliases weights_from_vector
#' @param membership_vector vector of strings containing group memberships
#' @param sep separator that delimits group memberships in membership_vector
#' @return sparse weight matrix
#' @export
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
  weights <- Matrix::Matrix(
    0,
    nrow=length(groups),
    ncol=length(membership_list),
    dimnames=list(as.character(groups), as.character(seq(length(membership_list))))
  )
  # fill in sparse weight matrix with 1s for each group membership
  weights[as.matrix(idx)] <- 1
  # return transposed sparse weight matrix
  return(weights)
}
