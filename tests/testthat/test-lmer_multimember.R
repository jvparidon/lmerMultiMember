test_that("return type is lmerMod", {
  df <- data.frame(
    x = seq(5),
    y = seq(5),
    memberships = c("a,b,c", "a,c", "a", "b", "b,a")
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer_multimember(y ~ x + (1|members), df, memb_mat=list(members=weights))
  expect_s4_class(m, "lmerMod")
})

test_that("lmer return type is lmerModMultiMember", {
  df <- data.frame(
    x = seq(5),
    y = seq(5),
    memberships = c("a,b,c", "a,c", "a", "b", "b,a")
  )
  weights <- weights_from_vector(df$memberships)
  m <- lmer_multimember(y ~ x + (1|members), df, memb_mat=list(members=weights))
  expect_s4_class(m, "lmerModMultiMember")
})

test_that("glmer return type is glmerModMultiMember", {
  df <- data.frame(
    x = runif(60, 0, 1),
    y = rbinom(60, 1, 0.6),
    memberships = rep(c("a,b,c", "a,c", "a", "b", "b,a", "b,c,a"), 10)
  )
  weights <- weights_from_vector(df$memberships)
  m <- glmer_multimember(y ~ x + (1|members), data=df, family=binomial, memberships=list(members=weights))
  expect_s4_class(m, "glmerModMultiMember")
})
