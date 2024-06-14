from dash_extensions.enrich import Dash,html,dcc,DashProxy,NoOutputTransform,CycleBreakerTransform,BlockingCallbackTransform,ServersideOutputTransform
import dash
import os
import sys
import pandas as pd
import overview as ov
import detail_graph as dg
import dash_bootstrap_components as dbc
import data_manager
import controller

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
app = DashProxy(__name__,transforms=[
            NoOutputTransform(),
            ],external_stylesheets=[dbc.themes.BOOTSTRAP],assets_ignore="lib/.*",use_pages=True)
main_page = dash.page_registry["pages.main_app"]
about_page = dash.page_registry["pages.about"]
dm = data_manager.DataManager(signatures,expressions,pathways)
ctrl = controller.Controller(dm)
df = pd.read_csv(signatures)
df["id"] = df["Cancer"]+"_"+df["Comparison"]
app.title = "THe_Biom"
detail_graph=dg.detail_graph()
mouse_up_event = {"event":"mouseup","props":["target","buttons","offsetX","offsetY","type"]}
mouse_down_event = {"event":"mousedown","props":["target","buttons","offsetX","offsetY","type"]}
mouse_out_event = {"event":"mouseout","props":["buttons","offsetX","offsetY","type"]}
mouse_move_event = {"event":"mousemove","props":["buttons","offsetX","offsetY","target","type"]}

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(
            html.H1("THe_Biom",id="title")
            ,width=11),
        dbc.Col(
                    [dcc.Link("About", href=about_page["relative_path"],style={"display":"auto"},id="about_link"),
                    dcc.Link("Main App", href=main_page["relative_path"],style={"display":"none"},id="main_link")],
                    width=1)
    ],style={"flexGrow":0,
             "flexShrink":1,
             }),
    dash.page_container
],fluid=True,style={"height":"100vh","display":"flex","flexDirection":"column"},id="main_container")
application = app.server
if __name__ == '__main__':
    # app.run(debug=False,host="172.17.11.246")
    app.run(debug=True)