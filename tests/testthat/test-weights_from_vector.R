test_that("returned data type is sparse matrix", {
  member_vec <-  c("a,b,c", "a,c", "a", "b", "b,a")
  expect_equal(is(weights_from_vector(member_vec), 'sparseMatrix'), TRUE)
})

test_that("returned sparse matrix has correct dimensions", {
  member_vec <-  c("a,b,c", "a,c", "a", "b", "b,a")
  expect_equal(dim(weights_from_vector(member_vec)), c(3, 5))
})

test_that("alternative string delimiters work", {
  member_vec <-  c("a;b;c", "a;c", "a", "b", "b;a")
  expect_equal(dim(weights_from_vector(member_vec, sep=";")), c(3, 5))
})
