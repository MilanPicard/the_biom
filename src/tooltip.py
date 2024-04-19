import dash_bootstrap_components as dbc
from dash import Dash, dcc
from dash import html


def create_tooltip(target):
    return dcc.Tooltip("",
                       id=f"{target.id}_tooltip",
                       show=False)
