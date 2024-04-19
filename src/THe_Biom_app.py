from dash import Dash,html,dcc
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
import tooltip

app = Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP],assets_ignore="lib/.*")
dm = data_manager.DataManager(sys.argv[1],sys.argv[2],sys.argv[3])
ctrl = controller.Controller(dm)
df = pd.read_csv(sys.argv[1])
df["id"] = df["Disease"]+"_"+df["Comparison"]
app.title = "THe_Biom"
detail_graph=dg.detail_graph()
app.layout = dbc.Container([
    dbc.Row([dbc.Col(
        menu.menu(dm.get_diseases(),dm.get_comparisons(),dm.get_disease_cmap(),dm.get_comparison_cmap(),dm.get_genes())
        ,width=2),
        dbc.Col([
            dbc.Row(dbc.Col(html.H1("THe_Biom",id="title"),width=12),style={"flexGrow":"0"}),
            dbc.Row([
                dbc.Col(ov.overview_graph(dm),width=6),
                dbc.Col([html.Span(tooltip.create_tooltip(detail_graph),style={"position":"relative","z-index":"4"}),dcc.Store(id="fake_graph_size"),detail_graph],width=6),
            ],style={
                "flexGrow":"1"
                ,"flexBasis":"40%"
                },
                className="g-0"),
            dbc.Row([
                # dcc.Tabs([
                # # dbc.Col(
                # dcc.Tab(label='Boxplot', children=[
                    dbp.detail_box_plot()
                #     # ,width=9)
                # ])
                #     ,
                # dcc.Tab(label='Heatmap', children=[

                # # dbc.Col(
                #     dbp.detail_heatmap()
                # ])
                    # ,width=3)
            ],style={
                "flexGrow":"1","flexBasis":"40%"
                },
                className="g-0")
        ],width=10,style={"display":"flex","flexDirection":"column","height":"100vh"}),
        html.Div(id="dummy_div"),
        html.Data(id="data_menu_selected"),
        html.Data(id="data_overview_selected"),
        html.Data(id="data_gene_detail_selected"),
        html.Data(id="data_gene_menu_selected"),
    ],style={"height":"100vh",})
],fluid=True)
if __name__ == '__main__':
    # app.run(debug=False,host="172.17.11.246")
    app.run(debug=True)