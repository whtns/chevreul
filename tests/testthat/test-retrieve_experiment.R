test_that("Assay retrieved", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
   mainExpName(chevreul_sce) <- "gene"
   retrieve_data <- retrieve_experiment(chevreul_sce, experiment = "gene")
    expect_equal(mainExpName(retrieve_data), "gene")
})
