test_that("possible to query experiment names", {
   chevreul_sce <- scuttle::mockSCE(ncells=200, ngenes=1000)
   mainExpName(chevreul_sce) <- "gene"
    expect_equal(query_experiment(chevreul_sce, "gene"), TRUE)
})
