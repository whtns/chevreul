test_that("Feature names obained", {
    names <- get_features(chevreul_sce)
    expect_equal(names, rownames(chevreul_sce))
})
