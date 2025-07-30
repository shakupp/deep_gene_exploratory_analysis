# Image de base avec R + Shiny Server
FROM rocker/shiny:latest

# Installer les dépendances système
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libhdf5-dev \
    libpng-dev \
    libglpk-dev \
    libudunits2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Copier les fichiers de l'application dans le conteneur
COPY . /app
WORKDIR /app

# Installer les packages CRAN
RUN R -e "install.packages(c( \
  'archive', 'colourpicker', 'corrplot', 'dplyr', 'DT', \
  'factoextra', 'fresh', 'ggplot2', 'ggrepel', 'glmnet', \
  'httr', 'jsonlite', 'lubridate', 'plotly', 'RColorBrewer', \
  'readr', 'reshape2', 'shiny', 'shinyalert', 'shinyBS', \
  'shinycssloaders', 'shinydashboard', 'shinydisconnect', \
  'shinyWidgets', 'survival' \
), repos = 'https://cloud.r-project.org')"

# Installer Bioconductor
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) \
  install.packages('BiocManager', repos='https://cloud.r-project.org'); \
  BiocManager::install(c( \
    'clusterProfiler', 'DESeq2', 'enrichplot', 'org.Mm.eg.db', \
    'pathview', 'ReactomePA' \
  ))"


RUN cp -r /app /srv/shiny-server/deep_gene_exploratory_analysis \
    && chown -R shiny:shiny /srv/shiny-server

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]