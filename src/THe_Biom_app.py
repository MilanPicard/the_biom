from dash_extensions.enrich import Dash,html,dcc
from dash_extensions import EventListener
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
if len(sys.argv)>=4:
    signatures = sys.argv[1]
    expressions = sys.argv[2]
    pathways = sys.argv[3]
elif ("THE_BIOM_SIGNATURES" in os.environ and "THE_BIOM_EXPRESSIONS" in os.environ and "THE_BIOM_PATHWAYS" in os.environ):
    signatures = os.environ["THE_BIOM_SIGNATURES"]
    expressions = os.environ["THE_BIOM_EXPRESSIONS"]
    pathways = os.environ["THE_BIOM_PATHWAYS"]
else:
    signatures =os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","signatures","THe_Biom_DEV_dataset.csv")
    expressions = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","activations","fake_data.csv")
    pathways = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","pathways","fake_pathways")
app = Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP],assets_ignore="lib/.*")
dm = data_manager.DataManager(signatures,expressions,pathways)
ctrl = controller.Controller(dm)
df = pd.read_csv(signatures)
df["id"] = df["Disease"]+"_"+df["Comparison"]
app.title = "THe_Biom"
detail_graph=dg.detail_graph()
mouse_up_event = {"event":"mouseup","props":["target","buttons","offsetX","offsetY","type"]}
mouse_down_event = {"event":"mousedown","props":["target","buttons","offsetX","offsetY","type"]}
mouse_out_event = {"event":"mouseout","props":["buttons","offsetX","offsetY","type"]}
mouse_move_event = {"event":"mousemove","props":["buttons","offsetX","offsetY","target","type"]}

app.layout = dbc.Container([
    dbc.Row([dbc.Col(
        menu.menu(dm.get_diseases(),dm.get_comparisons(),dm.get_disease_cmap(),dm.get_comparison_cmap(),dm.get_genes())
        ,width=2),
        EventListener([
            dbc.Row(dbc.Col(html.H1("THe_Biom",id="title"),width=12),style={"flexGrow":"0","height": "max-content"}),
            EventListener(
                [
                dbc.Col([
                    ov.overview_graph(dm),
                    html.Span(draggable="false",style={"width":"3%","height":"96%","cursor":"w-resize","userSelect":"none"}),
                         ],width=6,id="overview_col",style={"display":"flex","flexDirection":"row","borderWidth":"1px","borderStyle":"solid","borderColor":"black"}),
                dbc.Col([
                    html.Div([
                    html.Span(tooltip.create_tooltip(detail_graph),style={"position":"relative","zIndex":"4"}),
                         html.Span(id="detail_resize_span",draggable="false",style={"width":"3%","height":"100%","cursor":"w-resize","userSelect":"none"}),
                    dcc.Store(id="fake_graph_size",data={"AR":1}),
                         detail_graph,
                         dcc.Store(id="detail_graph_pos",data={}),
                         ],style={"display":"flex","flexDirection":"row","width":"100%","height":"96%"})
                         ],width=6,id="detail_col",style={"borderWidth":"1px","borderStyle":"solid","borderColor":"black"}),
            ],
            events=[mouse_move_event,mouse_up_event,mouse_down_event],useCapture=True,logging=False,id="move_in_ov",style={"flexGrow":"1","flexShrink":"1","flexBasis":"1%","height": "max-content","cursor":"n-resize"},className="g-0 row"),
            dbc.Row([
                # dcc.Tabs([
                # # dbc.Col(
                # dcc.Tab(label='Boxplot', children=[
                html.Div([#html.Span(style={"userSelect":"none","height":"1%"},id="height_resize_span"),
                          dbp.detail_box_plot()],style={"height":"100%","display":"flex","flexDirection":"column","justifyContent":"end"},id="second_row_div")
                #     # ,width=9)
                # ])
                #     ,
                # dcc.Tab(label='Heatmap', children=[

                # # dbc.Col(
                #     dbp.detail_heatmap()
                # ])
                    # ,width=3)
            ],style={
                "flexGrow":"1","flexBasis":"1%","flexShrink":"1","height": "max-content","borderWidth":"1px","borderStyle":"solid","borderColor":"black","cursor":"n-resize"
                },
                className="g-0",id="second_row")
        ],style={"display":"flex","flexDirection":"column","height":"100vh"},id="full_col",className="col-10 g-0",events=[mouse_down_event,mouse_up_event,mouse_move_event], logging=False),
        html.Div(id="dummy_div"),
        html.Data(id="data_menu_selected"),
        html.Data(id="data_overview_selected"),
        html.Data(id="data_gene_detail_selected"),
        html.Data(id="data_gene_menu_selected"),
        dcc.Store(id="resize_state",data={"width":{"is_resizing":False},"height":{"is_resizing":False}}),
    ],style={"height":"100vh","overflow":"hidden"},id="full_row")
],fluid=True)
application = app.server
if __name__ == '__main__':
    # app.run(debug=False,host="172.17.11.246")
    app.run(debug=True)