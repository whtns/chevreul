test_that("object metadata retrieved", {
    
    expect_type(get_object_metadata(chevreul_sce), "list")
})
