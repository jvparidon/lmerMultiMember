test_that("weights_from_vector() return type is sparse matrix", {
  a <- c("k,l,m", "k,m", "k", "l", "l,k")
  Wa <- weights_from_vector(a)
  expect_s4_class(Wa, "sparseMatrix")
})

test_that("weights_from_columns() return type is sparse matrix", {
  a <- cbind(
    c("k", "k", "l", NA, "m"),
    c("l", "m", "k", "m", "k"),
    c("m", NA, "m", NA, NA)
  )
  Wa <- weights_from_columns(a)
  expect_s4_class(Wa, "sparseMatrix")
})

test_that("interaction_weights() return type is sparse matrix", {
  a <- rep(c("k", "l", "l,k"), 2)
  b <- rep(c("m", "n"), 3)
  Wa <- weights_from_vector(a)
  Wb <- Matrix::fac2sparse(b)
  Wab <- interaction_weights(Wa, Wb)
  expect_s4_class(Wab, "sparseMatrix")
})

test_that("bradleyterry_from_vectors() return type is sparse matrix", {
  a <- c("k", "k", "l", "m", "m")
  b <- c("l", "m", "m", "k", "l")
  Wab <- bradleyterry_from_vectors(a, b)
  expect_s4_class(Wab, "sparseMatrix")
})

test_that("bradleyterry_from_sparse() return type is sparse matrix", {
  a <- Matrix::fac2sparse(c("k", "k", "l", "m", "m"))
  b <- Matrix::fac2sparse(c("l", "m", "m", "k", "l"))
  Wab <- bradleyterry_from_sparse(a, b)
  expect_s4_class(Wab, "sparseMatrix")
})

test_that("sparse matrix from weights_from_vector() has correct dimensions", {
  a <- c("k,l,m", "k,m", "k", "l", "l,k")
  Wa <- weights_from_vector(a)
  expect_equal(dim(Wa), c(3, 5))
})

test_that("sparse matrix from weights_from_columns() has correct dimensions", {
  a <- cbind(
    c("k", "k", "l", NA, "m"),
    c("l", "m", "k", "m", "k"),
    c("m", NA, "m", NA, NA)
  )
  Wa <- weights_from_columns(a)
  expect_equal(dim(Wa), c(3, 5))
})

test_that("sparse matrix from interaction_weights() has correct dimensions", {
  a <- rep(c("k", "l", "l,k"), 2)
  b <- rep(c("m", "n"), 3)
  Wa <- weights_from_vector(a)
  Wb <- Matrix::fac2sparse(b)
  Wab <- interaction_weights(Wa, Wb)
  expect_equal(dim(Wab), c(4, 6))
})

test_that("sparse matrix from bradleyterry_from_vectors() has correct dims", {
  a <- c("k", "k", "l", "m", "m")
  b <- c("l", "m", "m", "k", "l")
  Wab <- bradleyterry_from_vectors(a, b)
  expect_equal(dim(Wab), c(3, 5))
})

test_that("sparse matrix from bradleyterry_from_sparse() has correct dims", {
  a <- Matrix::fac2sparse(c("k", "k", "l", "m", "m"))
  b <- Matrix::fac2sparse(c("l", "m", "m", "k", "l"))
  Wab <- bradleyterry_from_sparse(a, b)
  expect_equal(dim(Wab), c(3, 5))
})

test_that("alternative string delimiters in weights_from_vector() work", {
  j <- c("a;b;c", "a;c", "a", "b", "b;a")
  Wj <- weights_from_vector(j, sep = ";")
  expect_equal(dim(Wj), c(3, 5))
})

test_that("returned sparse matrix from weights_from_vector() is correct", {
  Wreference <- as.matrix(t(data.frame(
    a = c(1, 1, 1, 0, 1),
    b = c(1, 0, 1, 0, 0),
    c = c(1, 1, 1, 1, 1),
    row.names = 1:5
  )))
  j <- c("a,b,c", "a,c", "b,a,c", "c", "c,a")
  Wj <- as.matrix(weights_from_vector(j))
  expect_equal(Wj, Wreference)
})

test_that("returned sparse matrix from weights_from_columns() is correct", {
  Wreference <- as.matrix(t(data.frame(
    a = c(1, 1, 1, 0, 1),
    b = c(1, 0, 1, 0, 0),
    c = c(1, 1, 1, 1, 1),
    row.names = 1:5
  )))
  j <- cbind(
    c("a", "a", "b", NA, "c"),
    c("b", "c", "a", "c", "a"),
    c("c", NA, "c", NA, NA)
  )
  Wj <- as.matrix(weights_from_columns(j))
  expect_equal(Wj, Wreference)
})

test_that("returned matrix from bradleyterry_from_vectors() is correct", {
  Wreference <- as.matrix(t(data.frame(
    k = c(1, 1, 0, -1, 0),
    l = c(-1, 0, 1, 0, -1),
    m = c(0, -1, -1, 1, 1),
    row.names = 1:5
  )))
  a <- c("k", "k", "l", "m", "m")
  b <- c("l", "m", "m", "k", "l")
  Wab <- as.matrix(bradleyterry_from_vectors(a, b))
  expect_equal(Wab, Wreference)
})

test_that("returned matrix from bradleyterry_from_sparse() is correct", {
  Wreference <- as.matrix(t(data.frame(
    k = c(1, 1, 0, -1, 0),
    l = c(-1, 0, 1, 0, -1),
    m = c(0, -1, -1, 1, 1),
    row.names = 1:5
  )))
  a <- Matrix::fac2sparse(c("k", "k", "l", "m", "m"))
  b <- Matrix::fac2sparse(c("l", "m", "m", "k", "l"))
  Wab <- as.matrix(bradleyterry_from_sparse(a, b))
  expect_equal(Wab, Wreference)
})
