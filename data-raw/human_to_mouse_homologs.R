# ## code to prepare `human_to_mouse_homologs` dataset goes here
#
# ## Not run:
# ## The code to prepare the .Rda file file from the marker file is:
# hs_genes <- rownames(chevreul_sce)
#
# biomaRt::listMarts(host="www.ensembl.org")
# human_mart <- biomaRt::useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
# mouse_mart <- biomaRt::useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
#
# human_to_mouse_homologs = biomaRt::getLDS(attributes = c("hgnc_symbol","entrezgene_id","ensembl_gene_id"),
#                                         filters = "hgnc_symbol", values = hs_genes,
#                                         # filters = "entrezgene_id", values = hs_entrez,
#                                         mart = human_mart,
#                                         attributesL = c("mgi_symbol","ensembl_gene_id","entrezgene_id"),
#                                         martL = mouse_mart)
#
# human_to_mouse_homologs <- human_to_mouse_homologs[c("HGNC.symbol", "MGI.symbol")]
#
# usethis::use_data(human_to_mouse_homologs, overwrite = TRUE)
