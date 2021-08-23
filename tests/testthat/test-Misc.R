source("test-common.R")

test_that("J_n, corresponding eigenvalues, an misc functions are okay",{
  
  expect_equal(J_n(5), J(5))
  
  expect_equal(eigenvals_J(3), eigenvals_J_n(3), ignore_attr = TRUE)
  
  
  expect_equal(mean_rice(140, 3), mean_rice_R(140, 3))
  expect_equal(mean_rice(70, 6), mean_rice_R(70, 6))
  expect_equal(mean_rice(20, 36), mean_rice_R(20, 36))

})


