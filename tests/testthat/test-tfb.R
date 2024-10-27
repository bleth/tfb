test_that("tfb-ols with wdim-att on sample df 1", {

  set.seed(2112)

  df <- tfb::tfb_sampledat1

  out <- tfb(df[,1:2],df[,3],df[,4],"ols","att")

  expect_equal(out$final$estimate,1.227715,tolerance = 1e-5)

})
