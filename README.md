# PEARL - Pathway and Gene Ontology Enrichment and Analysis Reporting Launcher

This project is an R Shiny application designed for the exploratory analysis of gene expression data. It provides interactive visualizations and analytical tools to facilitate the understanding of complex genomic datasets.

This application was developed as part of an R Shiny mock project for the BIMS (Bioinformatics) Master's program at the University of Rouen Normandy.

## Features

- **Interactive Visualizations**: Generate dynamic plots such as volcano plots, MA plots, and others to explore gene expression patterns.
- **Data Filtering**: Apply filters to focus on specific genes, conditions, or expression thresholds.
- **Statistical Analysis**: Perform differential expression analysis with customizable parameters.
- **User-Friendly Interface**: Navigate through the analysis with an intuitive and responsive UI.

## Installation

### Prerequisites

- R (version 4.0 or higher)
- RStudio (recommended)
- Git
- Conda (for environment management)
- Docker (for containerization)

### Clone the Repository

```bash
git clone https://github.com/shakupp/deep_gene_exploratory_analysis.git
cd deep_gene_exploratory_analysis

### Set Up Conda Environment (Optional)

Create a Conda environment using the provided environment.yml file:

conda env create -f environment.yml
conda activate deep_gene_env

### Run the Application Locally (without Docker)

R -e "shiny::runApp()"

### Dockerization

To containerize the application using Docker:
Build the Docker Image

docker build -t deep-gene-app .

### Run the Docker Container

docker run -p 3838:3838 deep-gene-app

Then access the application at:
http://localhost:3838/deep_gene_exploratory_analysis

## Citation

To cite this project, please refer to:

Sofr H., Martin Pena L., Maouloud L., Dauchel H. (2025). https://github.com/shakupp/deep_gene_exploratory_analysis

## License

This project is licensed under the MIT License. See the LICENSE file for details.