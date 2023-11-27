# shinyproxy
# FROM openanalytics/r-base
FROM openanalytics/r-ver:4.3.1

MAINTAINER Kevin Stachelek "kevin.stachelek@gmail.com"

RUN apt-get update && apt-get install --no-install-recommends -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libssl3 \
    libhdf5-dev \
    libboost-all-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    librsvg2-dev \
    cmake \
    build-essential \
    libglpk40 \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# system library dependency for the euler app
RUN apt-get update && apt-get install -y \
    libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("sf")'
RUN R -e "install.packages('igraph', dependencies = T, repos='https://cloud.r-project.org/')"
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install(version = "3.18")'
RUN R -e 'BiocManager::install(c("batchelor", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment"))'
RUN R -e 'BiocManager::install("pcaMethods")'
RUN R -e 'BiocManager::install("Biobase")'
RUN R -e 'BiocManager::install("tximport")'
RUN R -e 'BiocManager::install("annotables")'
RUN R -e 'BiocManager::install("genefilter")'
RUN R -e 'BiocManager::install("wiggleplotr")'
RUN R -e 'BiocManager::install("ensembldb")'
ARG GITHUB_PAT
ENV GITHUB_PAT=$GITHUB_PAT
RUN R -e 'install.packages("Seurat")'
RUN R -e 'install.packages("clustree", dependencies=TRUE)'
RUN R -e 'install.packages("ggpubr", dependencies=TRUE)'
RUN R -e 'remotes::install_github("velocyto-team/velocyto.R")'
RUN R -e 'remotes::install_github("cole-trapnell-lab/monocle3")'
#
# # RUN install2.r --error --deps TRUE unglue
# # RUN install2.r --error --deps TRUE rlang
# # RUN install2.r --error --deps TRUE fs
# # RUN install2.r --error --deps TRUE magrittr
# # RUN install2.r --error --deps TRUE future
# # RUN install2.r --error --deps TRUE DT
# # RUN install2.r --error --deps TRUE plotly
# # RUN install2.r --error --deps TRUE shinycssloaders
# # RUN install2.r --error --deps TRUE shinyFiles
# # RUN install2.r --error --deps TRUE janitor
# # RUN install2.r --error --deps TRUE shinyjqui
# # RUN install2.r --error --deps TRUE shinylogs
# # RUN install2.r --error --deps TRUE DBI
# # RUN install2.r --error --deps TRUE RSQLite
# # RUN install2.r --error --deps TRUE htmlwidgets
# # RUN install2.r --error --deps TRUE shinyWidgets
# # RUN install2.r --error --deps TRUE shinyjs
# # RUN install2.r --error --deps TRUE waiter
# # RUN install2.r --error --deps TRUE sessioninfo
# # RUN install2.r --error --deps TRUE httr
# # RUN install2.r --error --deps TRUE patchwork
# # RUN install2.r --error --deps TRUE textshaping
# # RUN install2.r --error --deps TRUE shiny
# # RUN install2.r --error --deps TRUE shinydashboard
#
#
# # COPY chevreul_*.tar.gz /app.tar.gz
# # RUN remotes::install_local('/app.tar.gz')
# # CMD R -e 'library(dockerfiler)'
#

# install dependencies of the euler app
RUN R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://cloud.r-project.org/')"
RUN R -e "install.packages('Rmpfr', repos='https://cloud.r-project.org/')"

# copy the app to the image
RUN mkdir /root/euler
COPY euler /root/euler

RUN R -e 'install.packages("R.utils", repos="https://cloud.r-project.org/")'
RUN R -e 'remotes::install_github("satijalab/seurat-wrappers")'

RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install matplotlib

RUN R -e 'install.packages("tidyverse")'
RUN R -e 'BiocManager::install("InteractiveComplexHeatmap")'

COPY Rprofile.site /usr/local/lib/R/etc/
RUN R -e 'remotes::install_github("whtns/chevreul")'
# EXPOSE 3838

# # install shinyproxy package with demo shiny application
# RUN apt-get update && apt-get install wget
# RUN wget https://github.com/openanalytics/shinyproxy-demo/raw/master/shinyproxy_0.0.1.tar.gz -O /root/shinyproxy_0.0# .1.tar.gz
# RUN R CMD INSTALL /root/shinyproxy_0.0.1.tar.gz
# RUN rm /root/shinyproxy_0.0.1.tar.gz
# # set host and port
# COPY application.yml /opt/shinyproxy/application.yml

COPY application.yml /opt/shinyproxy/application.yml

RUN mkdir /root/dockerapp
COPY dockerapp /root/dockerapp
CMD ["R", "-e", "shiny::runApp('/root/dockerapp')"]

FROM openanalytics/shinyproxy:3.0.0

COPY application.yml /opt/shinyproxy/application.yml

