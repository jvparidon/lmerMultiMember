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
