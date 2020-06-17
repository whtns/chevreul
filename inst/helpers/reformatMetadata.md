# Reformat Metadata

Arbitrary metadata can be appended to the data based on the results of exploratory data analysis. Metadata addition can be exectued by uploading a csv with rownames matching cell ids with new variables as columns.

seuratTools makes it simple to subset for a single batch or batch-integrated dataset. Subsetting can be accomplished either in a graphical setting by selection from a dimensionally reduced plot or by specifying a formatted file. Subsetiting of single or batch integrated data will trigger renewal of all relevant preprocessing steps including clustering, marker gene, and pathway enrichment as well as integration based on a 'batch' variable
