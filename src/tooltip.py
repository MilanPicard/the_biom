import dash_bootstrap_components as dbc
from dash_extensions.enrich import Dash, dcc
from dash_extensions.enrich import html


def create_tooltip(target):
    return dcc.Tooltip("",
                       id=f"{target.id}_tooltip",
                       show=False)
