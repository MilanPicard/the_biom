# THe Biom

Welcome to the development project repository of **THe Biom** (**T**CGA **H**EFS **Biom**arkers).

## Project Description

**THe_BIOM** is an interactive web application designed for the analysis and visualization of cancer biomarkers in transcriptomic data. It provides a comprehensive platform for researchers to explore gene signatures across different cancer types and disease stages. The application features an intuitive interface with multiple visualization tools, including graphs of signature relationships, detailed gene expression, and pathway enrichment visualization. Users can filter and analyze data by specific diseases, cancer stages, genes, and molecular pathways, enabling in-depth investigation of potential biomarkers. The platform also offers statistical analysis capabilities and various export options for data sharing and publication purposes.

## Prerequisites
Before you begin, ensure you have the following installed on your system:
- [Python 3.10](https://www.python.org/downloads/)
- [Git](https://git-scm.com/downloads/)

## Installation and Running the App

Below are platform-specific instructions to get the application running.

### For Windows Users

1.  **Open Command Prompt** (you can find it by searching for `cmd` in the Start Menu).

2.  **Clone the repository** and navigate into the project directory by running the following commands one by one:
    ```cmd
    git clone https://github.com/MilanPicard/the_biom.git
    cd the_biom
    ```

3.  **Create and activate a virtual environment**. This will create an isolated space for the project's dependencies. On Windows, the `py` launcher is the recommended way to ensure the correct Python version is used.
    ```cmd
    py -3.10 -m venv .venv
    .venv\Scripts\activate
    ```
    After running the second command, you should see `(.venv)` at the beginning of your command prompt line. If the first command fails, ensure you have Python 3.10 installed and that it's available in your system's PATH.

4.  **Install the required Python packages**:
    ```cmd
    pip install -r requirements.txt
    ```

5.  **Launch the application**:
    ```cmd
    start.bat
    ```

### For macOS and Linux Users

1.  **Open your terminal**.

2.  **Clone the repository** and navigate into the project directory by running the following commands:
    ```bash
    git clone https://github.com/MilanPicard/the_biom.git
    cd the_biom
    ```

3.  **Create and activate a virtual environment**:
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```
    After running the second command, you should see `(.venv)` at the beginning of your terminal prompt.
    *(Note: Depending on your system, you may need to use `python` instead of `python3`)*.

4.  **Install the required Python packages**:
    ```bash
    pip install -r requirements.txt
    ```

5.  **Launch the application**:
    ```bash
    bash start.sh
    ```

## Accessing the Application
Once the application is running, open your web browser and navigate to the address shown in the terminal (usually **http://127.0.0.1:8050/**).

**Note on Data Files:** The application uses hardcoded paths in the `start.sh` and `start.bat` files to load specific data sets from the `data/` directory. If you add new data, you will need to update these files to point to the new data files.


## Features
- Visualize gene signatures and their relationships across multiple cancer types and stages
- Interactive filtering by disease, stage, gene, and pathway
- Detailed gene expression and pathway enrichment views
- Export figures and data for publication or further analysis
- Download main datasets directly from the app
