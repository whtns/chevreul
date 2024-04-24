test_that("conversion works", {
  data("grch38_tx2gene")
  data("grch38")
    expect_contains(
        transcripts_to_genes(transcripts = c("ENST00000359842", "ENST00000470566", "ENST00000465764", "ENST00000619224")),
        c("RXRG")
    )
})
