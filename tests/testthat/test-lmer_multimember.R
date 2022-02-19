test_that("lmer return type is both lmerMod and lmerModMultiMember", {
  df <- data.frame(
    x = seq(5),
    y = seq(5),
    memberships = c("a,b,c", "a,c", "a", "b", "b,a")
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmerMultiMember::lmer(y ~ x + (1 | members),
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
  m <- lmerMultiMember::glmer(y ~ x + (1 | members),
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
  mm <- lmerMultiMember::lmer(Reaction ~ Days + (1 | Subject),
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
  mm <- lmerMultiMember::glmer(
    cbind(incidence, size - incidence) ~ period + (1 | herd),
    data = cbpp,
    family = binomial
  )
  expect_identical(fixef(l4), fixef(mm))
  expect_identical(ranef(l4), ranef(mm))
})

test_that("lmer with devFunOnly = TRUE works", {
  sleepstudy <- lme4::sleepstudy
  l4 <- lme4::lmer(Reaction ~ Days + (1 | Subject),
                   data = sleepstudy,
                   REML = FALSE,
                   devFunOnly = TRUE
  )
  mm <- lmerMultiMember::lmer(Reaction ~ Days + (1 | Subject),
                              data = sleepstudy,
                              REML = FALSE,
                              devFunOnly = TRUE
  )
  expect_identical(typeof(l4), typeof(mm))
})

test_that("glmer with devFunOnly = TRUE works", {
  cbpp <- lme4::cbpp
  l4 <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
                    data = cbpp,
                    family = binomial,
                    devFunOnly = TRUE
  )
  mm <- lmerMultiMember::glmer(
    cbind(incidence, size - incidence) ~ period + (1 | herd),
    data = cbpp,
    family = binomial,
    devFunOnly = TRUE
  )
  expect_identical(typeof(l4), typeof(mm))
})
