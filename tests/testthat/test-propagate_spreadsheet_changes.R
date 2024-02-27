test_that("Metadata updated", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  updated_table <- read_csv(system.file("extdata", "new_meta.csv", package="chevreul"))
  propagate_spreadsheet_changes(updated_table, chevreul_sce)

  expect_equal(colData(obj), updated_table)
})
