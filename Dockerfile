FROM openanalytics/r-base

MAINTAINER Kevin Stachelek "kstachelek@protonmail.com"

RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl1.0.0 \
    libhdf5-dev \
    libboost-all-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev

RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("sf")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install(c("batchelor", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment"))'
RUN R -e 'remotes::install_github("cole-trapnell-lab/monocle3")'
RUN R -e 'remotes::install_github("r-lib/remotes", ref = "6c8fdaa")'
RUN R -e 'install.packages("rlang")'
RUN R -e 'install.packages("fs")'
RUN R -e 'install.packages("rlang")'
RUN R -e 'install.packages("magrittr")'
RUN R -e 'BiocManager::install("pcaMethods")'
RUN R -e 'remotes::install_github("velocyto-team/velocyto.R")'
RUN R -e 'install.packages("future")'
RUN R -e 'install.packages("DT")'
RUN R -e 'install.packages("plotly")'
RUN R -e 'install.packages("shinycssloaders")'
RUN R -e 'install.packages("shinyFiles")'
RUN R -e 'install.packages("janitor")'
RUN R -e 'remotes::install_github("satijalab/seurat-wrappers")'
RUN R -e 'install.packages("unglue")'
RUN R -e 'install.packages("shinyjqui")'
RUN R -e 'BiocManager::install("Biobase")'
RUN R -e 'BiocManager::install("tximport")'
RUN R -e 'install.packages("clustree")'
RUN R -e 'install.packages("shinylogs")'
RUN R -e 'BiocManager::install("annotables")'
RUN R -e 'BiocManager::install("genefilter")'
RUN R -e 'remotes::install_github("immunogenomics/presto")'
RUN R -e 'install.packages("DBI")'
RUN R -e 'install.packages("RSQLite")'
RUN R -e 'install.packages("htmlwidgets")'
RUN R -e 'install.packages("shinyWidgets")'
RUN R -e 'install.packages("shinyjs")'
RUN R -e 'install.packages("waiter")'
RUN R -e 'install.packages("sessioninfo")'
RUN R -e 'remotes::install_github("mahmoudibrahim/genesorteR")'
RUN R -e 'install.packages("httr")'
RUN R -e 'install.packages("patchwork")'
RUN R -e 'remotes::install_github("jokergoo/ComplexHeatmap")'
RUN R -e 'BiocManager::install("wiggleplotr")'
RUN R -e 'BiocManager::install("ensembldb")'
RUN R -e 'install.packages("Seurat")'
RUN R -e 'install.packages("shiny")'
RUN R -e 'install.packages("shinydashboard")'
RUN R -e 'remotes::install_github("whtns/seuratTools")'

COPY seuratTools_*.tar.gz /app.tar.gz
RUN remotes::install_local('/app.tar.gz')
CMD R -e 'library(dockerfiler)'

# install shinyproxy package with demo shiny application
COPY shinyproxy_0.0.1.tar.gz /root/
RUN R CMD INSTALL /root/shinyproxy_0.0.1.tar.gz
RUN rm /root/shinyproxy_0.0.1.tar.gz

# set host and port
COPY Rprofile.site /usr/lib/R/etc/

EXPOSE 3838

CMD ["R", "-e", "shinyproxy::run_01_hello()"]
