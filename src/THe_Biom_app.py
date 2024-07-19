import dash
import os
import sys
import dash_bootstrap_components as dbc
import data_manager
if "THE_BIOM_MODE" in os.environ and os.environ["THE_BIOM_MODE"]=="PROD":
    from cel_app import celery_app
    background_callback_manager = dash.CeleryManager(celery_app)
elif "THE_BIOM_MODE" in os.environ and os.environ["THE_BIOM_MODE"]=="DEV_DISK":
    import diskcache
    cache = diskcache.Cache("./cache")
    background_callback_manager = dash.DiskcacheManager(cache)
else:
    background_callback_manager = None
app = dash.Dash(__name__,external_stylesheets=[dbc.themes.BOOTSTRAP],assets_ignore="lib/.*",use_pages=True,background_callback_manager=background_callback_manager)

about_page = dash.page_registry["pages.about"]
main_page = dash.page_registry["pages.main_app"]
app.title = "THe Biom"
mouse_up_event = {"event":"mouseup","props":["target","buttons","offsetX","offsetY","type"]}
mouse_down_event = {"event":"mousedown","props":["target","buttons","offsetX","offsetY","type"]}
mouse_out_event = {"event":"mouseout","props":["buttons","offsetX","offsetY","type"]}
mouse_move_event = {"event":"mousemove","props":["buttons","offsetX","offsetY","target","type"]}

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(
            dash.html.H1("THe Biom",id="title")
            ,width=11),
        dbc.Col(
                    [dash.dcc.Link("About", href=about_page["relative_path"],style={"display":"auto"},id="about_link"),
                    dash.dcc.Link("Main App", href=main_page["relative_path"],style={"display":"none"},id="main_link")],
                    width=1)
    ],style={"flexGrow":0,
             "flexShrink":1,
             }),
    dash.page_container
],fluid=True,style={"height":"100vh","display":"flex","flexDirection":"column"},id="main_container")
application = app.server
if __name__ == '__main__':
    import controller

    if len(sys.argv)>=4 and all([os.path.exists(i) for i in sys.argv[1:4]]):
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
    dm = data_manager.DataManager(signatures,expressions,pathways)
    controller.Controller.declare_callback(background_callback_manager!=None)
    ctrl = controller.Controller(dm)

    # app.run(debug=False,host="172.17.11.246")
    app.run(debug=True)