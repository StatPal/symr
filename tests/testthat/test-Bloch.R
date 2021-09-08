source("test-common.R")

test_that("Bloch equation is okay",{
  x <- c(56.27, 0.257, 0.56)
  TE <- c(0.5, 0.9); TR <- c(1, 5)
  expect_equal(Bloch_vec(x, TE, TR), Bloch_eqn_R(x, TE, TR), ignore_attr = TRUE)
  ## same values but may have different attributes (e.g., names and dimnames)
})


#test_that("Derivatives are okay",{
#  x <- c(56.27, 0.257, 0.56)
#  TE <- c(0.5, 0.9); TR <- c(1, 5)
#  
#  expect_equal(dee_v_ij_dee_W_ik(x, TE, TR, 1, 1), simple_dee_v_ij_dee_W_ik(x, TE, TR, 1, 1))
#  expect_equal(dee_v_ij_dee_W_ik(x, TE, TR, 1, 2), simple_dee_v_ij_dee_W_ik(x, TE, TR, 1, 2))
#  expect_equal(dee_v_ij_dee_W_ik(x, TE, TR, 1, 3), simple_dee_v_ij_dee_W_ik(x, TE, TR, 1, 3))
#}
#)

