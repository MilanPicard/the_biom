# from dash_extensions.enrich import html, register_page
from dash import html, register_page
from dash import dcc
import dash_bootstrap_components as dbc
register_page(__name__,"/about",title="About THe_BIOM")
def layout(**kwargs):
    # Define the tutorial steps as a list of dicts
    tutorial_steps = [
        {"src": "/assets/Mainvue.png", "header": "The main view of the app", "caption": "The filter menu on the left allows to focus on and visualize specific part of the data. You can select specific cancers or stages. They will be highlighted in the top left pannel, which shows the overview of the gene signatures as a graph. Each signature/stage is a node, while edges represent biomarkers overlap between two signatures."},
        {"src": "/assets/Overviewgraph.png", "header": "Overview graph", "caption": "The overview graph shows the available signatures (cancer/stage) as nodes, and the gene biomarkers in common as edges. Available cancers are UCEC(Uterine Corpus Endometrial Carcinoma), LUAD(Lung Adenocarcinoma), BRCA(Breast Invasive Carcinoma), LIHC(Liver Hepatocellular Carcinoma), KIRC(Kidney Renal Clear Cell Carcinoma), and HNSC(Head and Neck Squamous Cell Carcinoma). Clicking on a node selects a specific signature that will be shown in more details on the right pannel in Signature views. To select multiple signatures, use ctrl + click."},
        {"src": "/assets/Image2.png", "header": "Signature views", "caption": "Selecting signatures will open a detailed view in the top right pannel. Black dots represent biomarker genes and red triangles represent common pathways between signatures. Green nodes represent biomarkers that are present in multiple signatures and could correspond to important regulatory genes or predictive biomarkers. The husk around dots has a specific color depending on the cancer, and a specific striped pattern depending on the stage. Its position and form is arbitrary. You can click on legend items to select specific genes sets, and either copy them or send them to gProfiler for enrichment analysis."},
        {"src": "/assets/Image3.png", "header": "Gene selection", "caption": "By clicking on a gene (dots in the signature), it will open a detailed view of its expression across the different cancers and stages as boxplots. By selecting a gene in the filter menu on the left, it will, on top of showing boxplots expression values, open every signature in which this gene is considered a biomarker."},
        {"src": "/assets/Image_data.png", "header": "Add new signatures", "caption": "Users can add their own signatures to compare with the default signatures already present in the platform. By modifying the file ALLSIGNATURES_2024-07-08-14-21-18.csv in data/signatures and keeping the same format, they can add as many signatures (cancer/stage) as they want. Users can also remove existing signatures and only work with their signature data, THe Biom platform will automatically update the overview graph and signature views, and as such, be used as a generalizable visualisation tool."},
    ]
    return [
        dbc.Container([
            html.Div([
                html.H2("Development and exploitation of hybrid ensemble feature selection approaches: Application to the identification of cancer biomarkers in transcriptomic data", 
                    style={"fontSize": "1.5rem", "marginBottom": "1rem"}),
                html.P([
                    html.Strong("THe_BIOM"),
                    " is an interactive web application designed for the analysis and visualization of cancer biomarkers in transcriptomic data. It provides a comprehensive platform for researchers to explore gene signatures across different cancer types and disease stages. The application features an intuitive interface with multiple visualization tools, including graphs of signature relationships, detailed gene expression, and pathway enrichment visualization. Users can filter and analyze data by specific diseases, cancer stages, genes, and molecular pathways, enabling in-depth investigation of potential biomarkers. The platform also offers statistical analysis capabilities and various export options for data sharing and publication purposes."
                ], style={"fontSize": "1rem", "marginBottom": "1rem"}),
                html.Hr(),
                html.H2("Help", style={"fontSize": "1.5rem", "marginBottom": "1rem"}),
                # Tutorial steps as a vertical sequence
                *[
                    dbc.Card([
                        dbc.CardBody([
                            html.H4(step["header"], className="card-title", style={"marginBottom": "1rem", "fontWeight": "bold"}),
                            html.P(step["caption"], className="card-text", style={"fontSize": "1rem", "color": "#555", "marginBottom": "1.5rem"})
                        ]),
                        dbc.CardImg(src=step["src"], top=False, style={"maxWidth": "100%", "height": "auto", "margin": "0 auto", "display": "block", "borderRadius": "8px"}),
                    ], className="mb-5", style={"boxShadow": "0 2px 8px rgba(0,0,0,0.05)", "padding": "1rem", "background": "#fff"})
                    for step in tutorial_steps
                ],
                html.P("""Nature Communications: doi:10.1038/s41467-025-02400-0""", 
                    style={"fontSize": "1rem", "marginTop": "1rem"}),
                html.Footer(
                    dbc.Container([
                    dbc.Row(html.Hr()),
                    dbc.Row(html.P("""Acknowledgements"""), className="mb-3 mt-3"),
                    dbc.Row([
                        dbc.Col([
                            html.Img(src="/assets/Universite-Bordeaux-RVB-01_HD.jpg", style={"maxWidth": "100%", "height": "auto", "maxHeight": "100px", "display": "block", "margin": "0 auto"}),
                        ], width=2),
                        dbc.Col([
                            html.Img(src="/assets/UL.png", style={"maxWidth": "100%", "height": "auto", "maxHeight": "100px", "display": "block", "margin": "0 auto"}, alt="Logo université de Laval"),
                        ], width=2),
                        dbc.Col([
                            html.Img(src="/assets/adlab.png", style={"maxWidth": "100%", "height": "auto", "maxHeight": "100px", "display": "block", "margin": "0 auto"}, alt="Logo ADLab"),
                        ], width=2),
                        dbc.Col([
                            html.Img(src="/assets/CHUL.png", style={"maxWidth": "100%", "height": "auto", "maxHeight": "100px", "display": "block", "margin": "0 auto"}, alt="Logo CHUL"),
                        ], width=2),
                        dbc.Col([
                            html.Img(src="/assets/New_Logo_LaBRI.png", style={"maxWidth": "100%", "height": "auto", "maxHeight": "100px", "display": "block", "margin": "0 auto"}, alt="Logo LaBRI"),
                        ], width=2),
                    ], justify="center", align="center", className="mb-4 mt-2")
                ], fluid=True))
            ])
        ])
    ]
