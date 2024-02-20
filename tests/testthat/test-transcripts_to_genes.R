test_that("conversion works", {
  expect_error(
    transcripts_to_genes(transcripts = c("ENST00000359842", "ENST00000470566", "ENST00000465764", "ENST00000619224")),
    NA)
})
