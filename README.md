# THe Biom

Welcome to the development project repository of **THe Biom** (**T**CGA **H**EFS **Biom**arkers).

## Project Description

**THe_BIOM** is an interactive web application designed for the analysis and visualization of cancer biomarkers in transcriptomic data. It provides a comprehensive platform for researchers to explore gene signatures across different cancer types and disease stages. The application features an intuitive interface with multiple visualization tools, including graphs of signature relationships, detailed gene expression, and pathway enrichment visualization. Users can filter and analyze data by specific diseases, cancer stages, genes, and molecular pathways, enabling in-depth investigation of potential biomarkers. The platform also offers statistical analysis capabilities and various export options for data sharing and publication purposes.

## Requirements
- Python 3.10
- See `requirements.txt` for all dependencies

## Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/MilanPicard/the_biom.git
   cd the_biom
   ```
2. Install the necessary dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Running the App
To launch the application locally, run:
```bash
source .venv/bin/activate
bash start.sh
```

- By default, the app will use precalculated signatures located in the `data/` directory.
- The app will start a local web server. Open your browser and go to the address shown in the terminal (usually http://127.0.0.1:8050/).


## Features
- Visualize gene signatures and their relationships across multiple cancer types and stages
- Interactive filtering by disease, stage, gene, and pathway
- Detailed gene expression and pathway enrichment views
- Export figures and data for publication or further analysis
- Download main datasets directly from the app
