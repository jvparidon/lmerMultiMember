test_that("return type is lmerMod", {
  df <- data.frame(
    x = seq(60) + runif(60, 0, 10),
    y = seq(60) + rep(runif(6, 0, 10), 10),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer_multimember(y ~ x + (1|members), df, memberships=list(members=weights))
  expect_s4_class(m, "lmerMod")
})

test_that("return type is lmerModMultiMember", {
  df <- data.frame(
    x = seq(60) + runif(60, 0, 10),
    y = seq(60) + rep(runif(6, 0, 10), 10),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer_multimember(y ~ x + (1|members), df, memberships=list(members=weights))
  expect_s4_class(m, "lmerModMultiMember")
})
