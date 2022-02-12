source("test-common.R")

test_that("Bloch equation is okay", {
  x <- c(56.27, 0.257, 0.56)
  TE <- c(0.5, 0.9)
  TR <- c(1, 5)
  expect_equal(Bloch_vec(x, TE, TR), bloch(x, TE, TR), ignore_attr = TRUE)
  ## same values but may have different attributes (e.g., names and dimnames)
})


# test_that("Derivatives are okay",{
#  x <- c(56.27, 0.257, 0.56)
#  TE <- c(0.5, 0.9); TR <- c(1, 5)
#
#  expect_equal(dee_v_ij_dee_W_ik(x, TE, TR, 1, 1), simple_dee_v_ij_dee_W_ik(x, TE, TR, 1, 1))
#  expect_equal(dee_v_ij_dee_W_ik(x, TE, TR, 1, 2), simple_dee_v_ij_dee_W_ik(x, TE, TR, 1, 2))
#  expect_equal(dee_v_ij_dee_W_ik(x, TE, TR, 1, 3), simple_dee_v_ij_dee_W_ik(x, TE, TR, 1, 3))
# }
# )


test_that("Bloch equation as matrix is okay", {
  W <- rbind(c(50, 0.01, 0.003), c(36, 0.02, 0.04)) ## Two sample rows
  TE <- c(0.5, 0.9)
  TR <- c(1, 5)
  expect_equal(Bloch_vec_all(W, TE, TR), bloch.image(W, TE, TR), ignore_attr = TRUE)
  ## same values but may have different attributes (e.g., names and dimnames)
})


# test_that("Generated with noise is okay", {
#   W <- rbind(c(50, 0.01, 0.003), c(36, 0.02, 0.04)) ## Two sample rows
#   TE <- c(0.01, 0.03, 0.04, 0.01)
#   TR <- c(0.6, 0.6, 1, 0.8)
#   set.seed(1)
#   tmp1 = Generate.noisy.image(W, TE, TR, sig);
#   tmp2 = Generate.noisy.image_R(W, TE, TR, sig);
#   expect_equal(tmp1, tmp2, ignore_attr = TRUE)
#   ## same values but may have different attributes (e.g., names and dimnames)
# })