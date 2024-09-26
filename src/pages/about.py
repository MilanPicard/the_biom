# from dash_extensions.enrich import html, register_page
from dash import html, register_page
import dash_bootstrap_components as dbc
register_page(__name__,"/about",title="About THe_BIOM")
def layout(**kwargs):
    carousel = dbc.Carousel(
        items=[  
            {"key": "1","src":"/assets/start.png","caption":"","header":"This is the main view of the app"},
            {"key": "2","src":"/assets/start_overview.png","caption":"Each signature is a node.  The edges represent an overlap between two signatures","header":"The top left widget shows the overview of the gene signatures as a node-link diagram"},
            {"key": "3","src":"/assets/start_menu.png","caption":"","header":"The left menu is used to focus on specific part of the data"},
            {"key": "4","src":"/assets/start_diseases.png","caption":"","header":"You can select a subset of diseases you want to focus on"},
            {"key": "5","src":"/assets/start_comparisons.png","caption":"","header":"You can also focus on specific stage comparisons"},
            {"key": "6","src":"/assets/start_filter.png","caption":"","header":""},
            {"key": "7","src":"/assets/start_genes.png","caption":"","header":"You can select genes"},
            {"key": "8","src":"/assets/start_pathways.png","caption":"","header":"You can select pathways"},
            {"key": "9","src":"/assets/start_export.png","caption":"","header":"the export menu allows to save the representation or the underlying data of the different widgets"},
            {"key": "10","src":"/assets/mono_signature_view_with_cursor.png","caption":"","header":"Clicking on a signature on the overview widget show a detailled view of the signature in the mon signature widget."},
            {"key": "11","src":"/assets/multi_LIHC.png","caption":"Only the pathways involving genes belonging to multiple signatures are represented","header":"Selecting a non empty subset of diseases and comparisons shows a detailed view of the relevant signatures"},
            {"key": "12","src":"/assets/selecting_genes.png","caption":"You can select a gene either by using the left menu or by clicking on it in the signature views","header":"Selecting a gene focuses on the signatures this genes belongs to and shows expressions values(?) for this genes"},
            {"key": "13","src":"/assets/selecting_genes_filtering_diseases_stages.png","caption":"Only the diseases and comparisons selected on the menu are represented","header":"Selecting a gene focuses on the signatures this genes belongs to and shows expressions values(?) for this genes"},
            {"key": "14","src":"/assets/Selecting_pathways.png","caption":"","header":"Selecting a pathways in the menu focuses on the genes involved in this pathway"},
            {"key": "15","src":"/assets/stats.png","caption":"","header":"Overing a box shows detailed values as well as statistical significance between sets of conditions"},
        ],
        controls=True,indicators=True,class_name="carousel-fade",variant="dark"

    )
    return [
        html.Div([
            html.H2("Development and exploitation of hybrid ensemble feature selection approaches: Application to the identification of cancer biomarkers in transcriptomic data"),
            html.P("""Description of app : TODO"""),
            html.H2("Usage"),
            carousel,
            html.P("""submitted to ?? => DOI article si accepted"""),
            html.Footer(
                dbc.Container([
                dbc.Row(html.Hr()),
                dbc.Row(html.P("""Acknowledgements/logos tutelles Elsa""")),
                dbc.Row([
                    dbc.Col([
                    html.Img(src="/assets/Universite-Bordeaux-RVB-01_HD.jpg",style={"maxWidth":"100%"}),
                    ],width=3),
                    dbc.Col([
                    html.Img(src="/assets/UL.jpg",style={"maxWidth":"100%"},alt="Logo université de Laval"),
                    ],width=3),
                    dbc.Col([
                    html.Img(src="/assets/New_Logo_LaBRI.png",style={"maxWidth":"100%"},alt="Logo LaBRI?")
                    ],width=3)
                ])
            ]))
        ])
    ]
