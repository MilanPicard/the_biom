import sys
import os
import data_manager
# import controller
import pandas as pd
import detail_graph as dg
import detail_box_plot as dbp
import dash_bootstrap_components as dbc
import menu
from dash_extensions import EventListener
# from dash_extensions.enrich import Dash,html,dcc,register_page
from dash import Dash,html,dcc,register_page
import overview as ov
import tooltip
from dash import callback, Output, Input, State

# if len(sys.argv)>=4:
#     signatures = sys.argv[1]
#     expressions = sys.argv[2]
#     pathways = sys.argv[3]
# elif ("THE_BIOM_SIGNATURES" in os.environ and "THE_BIOM_EXPRESSIONS" in os.environ and "THE_BIOM_PATHWAYS" in os.environ):
#     signatures = os.environ["THE_BIOM_SIGNATURES"]
#     expressions = os.environ["THE_BIOM_EXPRESSIONS"]
#     pathways = os.environ["THE_BIOM_PATHWAYS"]
# else:
#     signatures =os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","signatures","THe_Biom_DEV_dataset.csv")
#     expressions = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","activations","fake_data.csv")
#     pathways = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),"data","pathways","fake_pathways.json")
register_page(__name__,"/",title="THe_BIOM")

# dm = data_manager.DataManager(signatures,expressions,pathways)
dm = data_manager.DataManager.get_instance()
# ctrl = controller.Controller(dm)
# df = pd.read_csv(signatures)
# df["id"] = df["Cancer"]+"_"+df["Comparison"]
mouse_up_event = {"event":"mouseup","props":["target","buttons","offsetX","offsetY","type"]}
click_event = {"event":"click","props":["target","buttons","offsetX","offsetY","type","isTrusted"]}
mouse_down_event = {"event":"mousedown","props":["target","buttons","offsetX","offsetY","type"]}
mouse_out_event = {"event":"mouseout","props":["buttons","offsetX","offsetY","type"]}
mouse_move_event = {"event":"mousemove","props":["buttons","offsetX","offsetY","target","type"]}
def layout():
    multi_graph=dg.detail_graph("detail_graph")
    mono_graph=dg.detail_graph("mono_graph")

    overview_graph = ov.overview_graph(dm)
    return dbc.Container([
        dbc.Row([
            dbc.Col(
                menu.menu(dm.get_diseases(),dm.get_comparisons(),dm.get_disease_cmap(),dm.get_comparison_cmap(),dm.get_genes(),dm.get_filters(),dm.get_pathways([],"Merge")),
                width=2,
                style={
                    "minWidth": '170px',
                    "height": "100%",
                    "overflow": "auto"
                }
            ),
            EventListener([
                EventListener(
                    [
                    dbc.Col([
                        html.Span(tooltip.create_tooltip(overview_graph),style={"position":"relative","zIndex":"4"}),
                        overview_graph,
                        dcc.Store(id="overview_graph_layout",data={}),
                        html.Span(draggable="false",style={"width":"3%","height":"96%","userSelect":"none"}),
                             ],width=6,id="overview_col",style={"display":"flex","flexDirection":"row","borderWidth":"1px","borderStyle":"solid","borderColor":"black","borderRadius":"5px"}),
                    dbc.Col([
                        html.Div([
                             html.Span(id="detail_resize_span",draggable="false",style={"width":"3%","height":"100%","userSelect":"none"},className="resize_span"),
                        dcc.Store(id="fake_graph_size",data={"AR":1}),
                        html.Div([
                        dbc.Tabs(id="detail_tabs",class_name="detail_tabs",active_tab="multi_signature_view",children=[
                            dbc.Tab(id="mono_tab",label="mono signature view",children=[
                        html.Span(tooltip.create_tooltip(mono_graph),style={"position":"relative","zIndex":"4"}),
                                
                                mono_graph,
                                                                                        html.Div(className="legend_canvas",children=[
                                                                                            html.Canvas(id="mono_canvas")
                                                                                            ])
                                                                                        ]),
                            dbc.Tab(tab_id="multi_signature_view",label="multi signature view",children=[
                        html.Span(tooltip.create_tooltip(multi_graph),style={"position":"relative","zIndex":"4"}),
                                
                                html.Div(
                                    style={"position": "relative", "width": "100%", "height": "100%"},
                                    children=[
                                        multi_graph,
                                        html.Div(
                                            id="multi_html_legend",
                                            className="legend_html_canvas",
                                            style={
                                                "position": "absolute",
                                                "top": "0px",
                                                "left": "0px",
                                                "padding": "2px 4px",
                                                "zIndex": 200,
                                                "minWidth": "0px",
                                                "maxWidth": "180px",
                                                "pointerEvents": "auto"
                                            }
                                        ),
                                        html.Div(
                                            id="selected_genes_legend",
                                            className="legend_html_canvas",
                                        ),
                                        dcc.Store(id="legend_tooltip_store"),
                                        dcc.Store(id="legend_hover_store"),
                                        dcc.Store(id="legend_active_store"),
                                    ]
                                ),
                                ]),
                            
                            ])],style={"flexGrow":1,"display":"flex","flexDirection":"column"})
                             ,
                             dcc.Store(id="detail_graph_pos",data={}),
                             dcc.Store(id="mono_graph_pos",data={}),
                             dcc.Store(id="mono_legend",data={}),
                             dcc.Store(id="multi_legend",data={}),
                             dcc.Store(id="multi_export",data={}),
                             dcc.Store(id="mono_export",data={}),
                             dcc.Store(id="clipboard_text_store",data=""),
                             dcc.Store(id="gprofiler_url_store",data=""),
                             dcc.Store(id="copy_toast_store",data={}),
                             dbc.Toast(
                                 id="copy_toast",
                                 header="Copy to Clipboard",
                                 is_open=False,
                                 duration=3000,
                                 icon="success",
                                 style={"position": "fixed", "bottom": 20, "right": 20, "zIndex": 9999},
                                 children=""
                             ),
                             ],style={"display":"flex","flexDirection":"row","width":"100%","height":"96%"})
                             ],width=6,id="detail_col",style={"borderWidth":"1px","borderStyle":"solid","borderColor":"black","borderRadius":"5px"}),
                ],
                events=[click_event,mouse_move_event,mouse_up_event,mouse_down_event],useCapture=True,logging=False,id="move_in_ov",style={"flex": "0.95 1 0","minHeight": "200px","display": "flex","flexDirection": "row"},className="g-0 row"),
                dbc.Row([
                    html.Div(
                        dbp.detail_box_plot(),
                        style={
                            "height": "100%",
                            "display": "flex",
                            "flexDirection": "column",
                            "overflowY": "auto",
                            "minHeight": "0",
                            "flex": "1 1 auto"
                        },
                        id="second_row_div"
                    )
                ], style={
                    "flex": "1.05 1 0",
                    "minHeight": "120px",
                    "maxHeight": "100%",
                    "display": "flex",
                    "flexDirection": "column",
                    "borderWidth": "1px",
                    "borderStyle": "solid",
                    "borderColor": "black",
                    "borderRadius": "5px",
                    "overflow": "hidden"
                },
                className="g-0",
                id="second_row")
            ],
            style={
                "display": "flex",
                "flexDirection": "column",
                "height": "100%",
                "maxWidth": 'calc(100vw - 170px)',
                "overflow": "hidden"
            },
            className="col-10 g-0",
            events=[mouse_down_event,mouse_up_event,mouse_move_event], logging=False),
            html.Div(id="dummy_div",style={"flexBasis":"0px","flexGrow":0}),
            html.Data(id="data_menu_selected",style={"flexBasis":"0px","flexGrow":0}),
            html.Data(id="data_overview_selected",style={"flexBasis":"0px","flexGrow":0}),
            html.Data(id="data_gene_detail_selected",style={"flexBasis":"0px","flexGrow":0}),
            html.Data(id="data_gene_menu_selected",style={"flexBasis":"0px","flexGrow":0}),
            dcc.Store(id="resize_state",data={"width":{"is_resizing":False},"height":{"is_resizing":False}}),
        ],
        style={
            "height": "100%",
            "margin": "0",
            "padding": "0",
            "display": "flex",
            "flexDirection": "column",
            "overflow": "hidden"
        },
        className="g-0")
    ],
    fluid=True,
    style={
        "height": "100%",
        "padding": "0",
        "margin": "0",
        "maxWidth": "100vw"
    })
# ,style={
#         "flexGrow":1,
#         "flexShrink":0,
#         "flewBasis":"80vh",
#         #"height":"100vh",
#              "overflow":"hidden"},id="full_row")

@callback(
    Output("copy_toast", "is_open"),
    Output("copy_toast", "children"),
    Input("copy_toast_store", "data"),
    State("copy_toast", "is_open"),
    prevent_initial_call=True
)
def update_copy_toast(toast_data, is_open):
    if toast_data and toast_data.get("is_open"):
        return True, toast_data.get("message", "Successfully copied genes")
    return False, ""
