test_that("lmer return type is both lmerMod and lmerModMultiMember", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer(y ~ x + (1 | members),
    data = df,
    memberships = list(members = weights)
  )
  expect_s4_class(m, "lmerMod")
  expect_s4_class(m, "lmerModMultiMember")
})

test_that("glmer return type is both glmerMod and glmerModMultiMember", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- glmer(y ~ x + (1 | members),
    data = df,
    family = binomial,
    memberships = list(members = weights)
  )
  expect_s4_class(m, "glmerMod")
  expect_s4_class(m, "glmerModMultiMember")
})


test_that("lmer pass through to lme4 works", {
  sleepstudy <- lme4::sleepstudy
  l4 <- lme4::lmer(Reaction ~ Days + (1 | Subject),
    data = sleepstudy,
    REML = FALSE
  )
  mm <- lmer(Reaction ~ Days + (1 | Subject),
    data = sleepstudy,
    REML = FALSE
  )
  expect_identical(fixef(l4), fixef(mm))
  expect_identical(ranef(l4), ranef(mm))
  expect(!lme4::isREML(mm), "REML passthrough failed")
})

test_that("glmer pass through to lme4 works", {
  cbpp <- lme4::cbpp
  l4 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
    data = cbpp,
    family = binomial
  )
  mm <- glmer(
    cbind(incidence, size - incidence) ~ period + (1 | herd),
    data = cbpp,
    family = binomial
  )
  expect_identical(fixef(l4), fixef(mm))
  expect_identical(ranef(l4), ranef(mm))
})

test_that("lmer with devFunOnly = TRUE works", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer(y ~ x + (1 | members),
    data = df,
    memberships = list(members = weights),
    devFunOnly = TRUE
  )
  expect_identical(typeof(m), "closure")
})

test_that("glmer with devFunOnly = TRUE works", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- glmer(y ~ x + (1 | members),
    data = df,
    family = binomial,
    memberships = list(members = weights),
    devFunOnly = TRUE
  )
  expect_identical(typeof(m), "closure")
})

test_that("glmer with devFunOnly = TRUE & nAGQ = 0 works", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- glmer(y ~ x + (1 | members),
    data = df,
    family = binomial,
    memberships = list(members = weights),
    nAGQ = 0,
    devFunOnly = TRUE
  )
  expect_identical(typeof(m), "closure")
})

test_that("glmer with start = list(theta = .8) works", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  expect_silent(glmer(y ~ x + (1 | members),
    data = df,
    family = binomial,
    memberships = list(members = weights),
    start = list(theta = .8)
  ))
})

test_that("glmer with start = .8 works", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  expect_silent(glmer(y ~ x + (1 | members),
    data = df,
    family = binomial,
    memberships = list(members = weights),
    start = .8
  ))
})

test_that("glmer with invalid start parameter throws error", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  expect_error(glmer(y ~ x + (1 | members),
    data = df,
    family = binomial,
    memberships = list(members = weights),
    start = list(heron = .8)
  ))
})

test_that("glmer with fixef start parameter & nAGQ = 0 throws error", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  expect_error(glmer(y ~ x + (1 | members),
    data = df,
    family = binomial,
    memberships = list(members = weights),
    nAGQ = 0,
    start = list(fixef = c(.2, .8))
  ))
})

test_that("summary.merModMultiMember wih incorrect type throws error", {
  expect_error(summary.merModMultiMember(0))

  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  member_matrix <- weights_from_vector(df$memberships)
  m <- lmer(y ~ x + (1 | members),
                             data = df,
                             memberships = list(members = member_matrix)
  )
  class(m) <- c("merModMultiMember")
  expect_error(summary(m))
})

test_that("summary.merModMultiMember has the right attributes and prints", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  member_matrix <- weights_from_vector(df$memberships)
  # call with lmer
  m <- lmer(y ~ x + (1 | members),
    data = df,
    memberships = list(members = member_matrix)
  )
  expect_equal(member_matrix, summary(m)$memberships$members)
  expect_output(print(summary(m)))
  # call with glmer
  m <- glmer(y ~ x + (1 | members),
                              data = df,
                              family = binomial,
                              memberships = list(members = member_matrix)
  )
  expect_equal(member_matrix, summary(m)$memberships$members)
  expect_output(print(summary(m)))
})

test_that("plot_membership_hist returns a plot for each object type", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  member_matrix <- weights_from_vector(df$memberships)
  expect_s3_class(plot_membership_hist(member_matrix),
                  "histogram")
  # call with lmer
  m <- lmer(y ~ x + (1 | members),
                             data = df,
                             memberships = list(members = member_matrix)
  )
  expect_s3_class(plot_membership_hist(m, "members"),
                  "histogram")
  expect_s3_class(plot_membership_hist(summary(m), "members"),
                  "histogram")
})

test_that("plot_membership_matrix returns a plot for each object type", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  member_matrix <- weights_from_vector(df$memberships)
  expect_s3_class(plot_membership_matrix(member_matrix),
                  "trellis")
  # call with lmer
  m <- lmer(y ~ x + (1 | members),
                             data = df,
                             memberships = list(members = member_matrix)
  )
  expect_s3_class(plot_membership_matrix(m, "members"),
                  "trellis")
  expect_s3_class(plot_membership_matrix(summary(m), "members"),
                  "trellis")
})
