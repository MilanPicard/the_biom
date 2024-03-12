from dash import Dash,html
import os
import sys
import pandas as pd
import overview as ov
import detail_graph as dg
import interactions
import detail_box_plot as dbp
import dash_bootstrap_components as dbc
import data_manager
import controller
import menu

app = Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP])
dm = data_manager.DataManager(sys.argv[1],sys.argv[2],sys.argv[3])
ctrl = controller.Controller(dm)
df = pd.read_csv(sys.argv[1])
df["id"] = df["Disease"]+"_"+df["Comparison"]
app.title = "THe_Biom"
app.layout = dbc.Container([
    dbc.Row([dbc.Col(
        menu.menu(dm.get_diseases(),dm.get_stages(),dm.get_disease_cmap(),dm.get_stage_cmap())
        ,width=2),
        dbc.Col([
            dbc.Row(dbc.Col(html.H1("THe_Biom",id="title"),width=12),style={"flex-grow":"0"}),
            dbc.Row([
                dbc.Col(ov.overview_graph(dm),width=6),
                dbc.Col(dg.detail_graph(),width=6),
            ],style={
                "flex-grow":"1"
                ,"flex-basis":"40%"
                }),
            dbc.Row([
                dbc.Col(dbp.detail_box_plot(),width=9),
                dbc.Col(dbp.detail_heatmap(),width=3)
            ],style={
                "flex-grow":"1","flex-basis":"40%"
                })
        ],width=10,style={"display":"flex","flex-direction":"column","height":"100vh"}),
        html.Div(id="dummy_div")
    ],style={"height":"100vh",})
],fluid=True)
if __name__ == '__main__':
    app.run(debug=False,host="172.17.11.246")