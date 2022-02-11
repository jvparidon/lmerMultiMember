test_that("return type is lmerMod", {
  df <- data.frame(
    x = seq(5),
    y = seq(5),
    memberships = c("a,b,c", "a,c", "a", "b", "b,a")
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer(y ~ x + (1|members), df, memberships=list(members=weights))
  expect_s4_class(m, "lmerMod")
})

test_that("return type is lmerModMultiMember", {
  df <- data.frame(
    x = seq(5),
    y = seq(5),
    memberships = c("a,b,c", "a,c", "a", "b", "b,a")
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer(y ~ x + (1|members), df, memberships=list(members=weights))
  expect_s4_class(m, "lmerModMultiMember")
})


test_that("pass through to lme4 works", {
  sleepstudy <- lme4::sleepstudy
  l4 <- lme4::lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
  mm <- lmerMultiMember::lmer(Reaction ~ Days + (1|Subject), sleepstudy, REML=FALSE)
  expect_identical(fixef(l4), fixef(mm))
  expect_identical(ranef(l4), ranef(mm))
  expect(!lme4::isREML(mm), "REML passthrough failed")
})
