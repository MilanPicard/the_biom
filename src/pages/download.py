from dash import html, register_page
import dash_bootstrap_components as dbc
register_page(__name__, "/download", title="Download Data")
def layout(**kwargs):
    return [
        dbc.Container([
            html.H2("Download Data", style={"fontSize": "1.25rem", "marginBottom": "0.5rem", "marginTop": "1rem"}),
            html.P("You can download the main datasets used in this application for your own analysis:", style={"fontSize": "1rem", "marginBottom": "0.5rem"}),
            dbc.ListGroup([
                dbc.ListGroupItem([
                    html.Span("Gene Expression Data (Activations): "),
                    html.A("Download CSV (12 MB)", href="/download/data/activations/ALLEXPRESSION_log2TPMplus1_2024-07-09-16-19-05.csv", download="ALLEXPRESSION_log2TPMplus1.csv", target="_blank", style={"marginLeft": "0.5rem"})
                ]),
                dbc.ListGroupItem([
                    html.Span("Pathways Data: "),
                    html.A("Download JSON (3.3 MB)", href="/download/data/pathways/pathways.json", download="pathways.json", target="_blank", style={"marginLeft": "0.5rem"})
                ]),
                dbc.ListGroupItem([
                    html.Span("Signatures Data: "),
                    html.A("Download CSV (22 KB)", href="/download/data/signatures/ALLSIGNATURES_2024-07-08-14-21-18.csv", download="ALLSIGNATURES.csv", target="_blank", style={"marginLeft": "0.5rem"})
                ]),
            ], flush=True, style={"marginBottom": "1.5rem"}),
            html.Div([
                html.H3("Source Code", style={"fontSize": "1.25rem", "marginTop": "2rem", "marginBottom": "0.5rem"}),
                html.P([
                    "You can find the source code and contribute to the project on ",
                    html.A("MilanPicard/the_biom", href="https://github.com/MilanPicard/the_biom", target="_blank", style={"color": "#0366d6", "textDecoration": "underline"}),
                    "."
                ], style={"fontSize": "1rem", "marginBottom": "1rem"})
            ]),
            html.Div([
                html.H3("References", style={"fontSize": "1.25rem", "marginTop": "2rem", "marginBottom": "0.5rem"}),
                html.Ul([
                    html.Li([
                        "THe Biom: a database of cancers novel transcriptomic biomarkers identified by robust feature selection"
                    ]),
                    html.Li([
                        "Claude, E., Leclercq, M., Thébault, P., Droit, A., & Uricaru, R. (2024). Optimizing hybrid ensemble feature selection strategies for transcriptomic biomarker discovery in complex diseases. NAR Genomics and Bioinformatics, 6(3), lqae079."
                    ])
                ], style={"fontSize": "1rem", "marginBottom": "1rem"})
            ]),
        ])
    ] 