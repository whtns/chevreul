test_that("returns a character vector", {
    
    expect_type(
        metadata_from_object(chevreul_sce),
        "character"
    )
})
