source("test-common.R")

## Implement and check the function with mask


test_that("Performance measures are okay",{
    W <- rbind(c(40, 0.01, 0.003), c(36, 0.02, 0.04)) ## Two sample rows of parameters, W
    test <- rbind(c(56, 52, 57, 51), c(39, 37, 33, 34.4))

    ## Test design parameters
    TE <- c(0.01, 0.03, 0.04, 0.01)
    TR <- c(0.6, 0.6, 1, 0.8)
    sig <- c(1.2, 1.1, 1.4, 1.2)
    mask <- c(0, 0)   # make the other test function with mask

    tmp11 = performance(W, test, TE, TR, sig, mask, 1, 1, 1)
    tmp12 = performance(W, test, TE, TR, sig, mask, 3, 1, 1)
    tmp13 = performance(W, test, TE, TR, sig, mask, 1, 2, 1)
    tmp14 = performance(W, test, TE, TR, sig, mask, 3, 2, 1)

    tmp21 = Performance_test(W, test, TE, TR, sig, 1, 1, 1)
    tmp22 = Performance_test(W, test, TE, TR, sig, 3, 1, 1)
    tmp23 = Performance_test(W, test, TE, TR, sig, 1, 2, 1)
    tmp24 = Performance_test(W, test, TE, TR, sig, 3, 2, 1)


    expect_equal(tmp11, tmp21)
    expect_equal(tmp12, tmp22)
    expect_equal(tmp13, tmp23)
    expect_equal(tmp14, tmp24)
}
)



test_that("Performance measure absolute value with prediction",{
    W <- rbind(c(40, 0.01, 0.003), c(36, 0.02, 0.04)) ## Two sample rows of parameters, W
    test <- rbind(c(56, 52, 57, 51), c(39, 37, 33, 34.4))

    ## Test design parameters
    TE <- c(0.01, 0.03, 0.04, 0.01)
    TR <- c(0.6, 0.6, 1, 0.8)
    sig <- c(1.2, 1.1, 1.4, 1.2)
    mask <- c(0, 0)   # make the other test function with mask

    img = bloch.image(W, TE, TR);

    tmp11 = performance(W, test, TE, TR, sig, mask, 1, 1, 0)
    tmp13 = performance(W, test, TE, TR, sig, mask, 1, 2, 0)

    expect_equal(tmp11, colMeans(abs(img - test)))
    expect_equal(tmp13, sqrt(colMeans((img - test)^2)))
}
)
