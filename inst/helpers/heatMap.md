# Heatmap

Heatmap of expression using 'Seurat::DoHeatMap'. 

Expression values for each cell are normalized by that cell's total expression then multiplied by 10,000 and natural-log transformed before plotting.

By default, 50 most highly variable genes are displayed. However, an arbitrary lists of genes can be plotted for comparison. The genes/transcripts are displayed in the order they are listed.
