test_that("return type is sparse matrix", {
  member_vec <-  c("a,b,c", "a,c", "a", "b", "b,a")
  expect_s4_class(weights_from_vector(member_vec), 'sparseMatrix')

  member_cols <- cbind(
    c("a", "a", "b", NA, "c"),
    c("b", "c", "a", "c", "a"),
    c("c", NA, "c", NA, NA)
  )
  expect_s4_class(weights_from_columns(member_cols), 'sparseMatrix')
})

test_that("returned sparse matrix has correct dimensions", {
  member_vec <-  c("a,b,c", "a,c", "a", "b", "b,a")
  expect_equal(dim(weights_from_vector(member_vec)), c(3, 5))

  member_cols <- cbind(
    c("a", "a", "b", NA, "c"),
    c("b", "c", "a", "c", "a"),
    c("c", NA, "c", NA, NA)
  )
  expect_equal(dim(weights_from_columns(member_cols)), c(3, 5))
})

test_that("alternative string delimiters work", {
  member_vec <-  c("a;b;c", "a;c", "a", "b", "b;a")
  expect_equal(dim(weights_from_vector(member_vec, sep=";")), c(3, 5))
})

test_that("returned sparse matrices are correct", {
  reference <- as.matrix(t(data.frame(a = c(1, 1, 1, 0, 1),
                                      b = c(1, 0, 1, 0, 0),
                                      c = c(1, 1, 1, 1, 1),
                                      row.names = 1:5)))

  test_vec <- as.matrix(weights_from_vector(member_vec <- c(
    "a,b,c",
    "a,c",
    "b,a,c",
    "c",
    "c,a"
  )))
  expect_equal(test_vec, reference)

  test_cols <- as.matrix(weights_from_columns(cbind(
    c("a", "a", "b", NA, "c"),
    c("b", "c", "a", "c", "a"),
    c("c", NA, "c", NA, NA)
  )))
  expect_equal(test_cols, reference)
})
