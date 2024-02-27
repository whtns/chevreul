test_that("Project renamed", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  obj<-rename_object(chevreul_sce, "new_name")
  expect_contains(metadata(obj)["project.name"], "new_name")
})
