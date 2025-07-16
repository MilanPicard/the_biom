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
        ])
    ] 