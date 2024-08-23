test_that("Project renamed", {
    
    obj <- rename_object(chevreul_sce, "new_name")
    expect_contains(metadata(obj)["project.name"], "new_name")
})
