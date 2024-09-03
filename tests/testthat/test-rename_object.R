test_that("Project renamed", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
    obj <- rename_object(chevreul_sce, "new_name")
    expect_contains(metadata(obj)["project.name"], "new_name")
})
