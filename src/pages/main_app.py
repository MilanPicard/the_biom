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
    return dbc.Container([ dbc.Row(
     [
        dbc.Col(
        menu.menu(dm.get_diseases(),dm.get_comparisons(),dm.get_disease_cmap(),dm.get_comparison_cmap(),dm.get_genes(),dm.get_filters(),dm.get_pathways([],"Merge"))
        ,width=2,style={"minWidth":'170px',"height":"100%"}),
        EventListener([
            # dbc.Row([
            #     dbc.Col(html.H1("THe_Biom",id="title"),width=11),
            #     dbc.Col(
            #         # dcc.Link("About", href=about_page["relative_path"]),
            #         width=1)
            #          ],style={"flexGrow":"0","height": "max-content"}),
            EventListener(
                [
                dbc.Col([
                    html.Span(tooltip.create_tooltip(overview_graph),style={"position":"relative","zIndex":"4"}),
                    overview_graph,
                    dcc.Store(id="overview_graph_layout",data={}),
                    html.Span(draggable="false",style={"width":"3%","height":"96%","cursor":"w-resize","userSelect":"none"}),
                         ],width=6,id="overview_col",style={"display":"flex","flexDirection":"row","borderWidth":"1px","borderStyle":"solid","borderColor":"black","borderRadius":"5px"}),
                dbc.Col([
                    html.Div([
                         html.Span(id="detail_resize_span",draggable="false",style={"width":"3%","height":"100%","cursor":"w-resize","userSelect":"none"},className="resize_span"),
                    dcc.Store(id="fake_graph_size",data={"AR":1}),
                    html.Div([
                    dbc.Tabs(id="detail_tabs",class_name="detail_tabs",children=[
                        dbc.Tab(id="mono_tab",label="mono signature view",children=[
                    html.Span(tooltip.create_tooltip(mono_graph),style={"position":"relative","zIndex":"4"}),
                            
                            mono_graph,
                                                                                    html.Div(className="legend_canvas",children=[
                                                                                        html.Canvas(id="mono_canvas")
                                                                                        ])
                                                                                    ]),
                        dbc.Tab(label="multi signature view",children=[
                    html.Span(tooltip.create_tooltip(multi_graph),style={"position":"relative","zIndex":"4"}),
                            
                            multi_graph,html.Div(className="legend_canvas",children=[html.Canvas(id="multi_canvas")])]),
                        
                        ])],style={"flexGrow":1,"display":"flex","flexDirection":"column"})
                         ,
                         dcc.Store(id="detail_graph_pos",data={}),
                         dcc.Store(id="mono_graph_pos",data={}),
                         dcc.Store(id="mono_legend",data={}),
                         dcc.Store(id="multi_legend",data={}),
                         dcc.Store(id="multi_export",data={}),
                         dcc.Store(id="mono_export",data={}),
                         ],style={"display":"flex","flexDirection":"row","width":"100%","height":"96%"})
                         ],width=6,id="detail_col",style={"borderWidth":"1px","borderStyle":"solid","borderColor":"black","borderRadius":"5px"}),
            ],
            events=[click_event,mouse_move_event,mouse_up_event,mouse_down_event],useCapture=True,logging=False,id="move_in_ov",style={"flexGrow":"1","flexShrink":"1","flexBasis":"1%","height": "max-content","cursor":"n-resize"},className="g-0 row"),
            dbc.Row([
                # dcc.Tabs([
                # # dbc.Col(
                # dcc.Tab(label='Boxplot', children=[
                html.Div([#html.Span(style={"userSelect":"none","height":"1%"},id="height_resize_span")
                    ]+ dbp.detail_box_plot(),style={"height":"100%","display":"flex","flexDirection":"column","justifyContent":"end"},id="second_row_div")
                #     # ,width=9)
                # ])
                #     ,
                # dcc.Tab(label='Heatmap', children=[

                # # dbc.Col(
                #     dbp.detail_heatmap()
                # ])
                    # ,width=3)
            ],style={
                "flexGrow":"1","flexBasis":"1%","flexShrink":"1","height": "max-content","borderWidth":"1px","borderStyle":"solid","borderColor":"black","cursor":"n-resize","borderRadius":"5px","overflowY":"auto"
                },
                className="g-0",id="second_row")
        ],style={"display":"flex","flexDirection":"column","height":"100%","maxWidth":'calc(100% - 170px)'
                #  ,"height":"100vh"
                 },id="full_col",className="col-10 g-0",events=[mouse_down_event,mouse_up_event,mouse_move_event], logging=False),
        html.Div(id="dummy_div",style={"flexBasis":"0px","flexGrow":0}),
        html.Data(id="data_menu_selected",style={"flexBasis":"0px","flexGrow":0}),
        html.Data(id="data_overview_selected",style={"flexBasis":"0px","flexGrow":0}),
        html.Data(id="data_gene_detail_selected",style={"flexBasis":"0px","flexGrow":0}),
        html.Data(id="data_gene_menu_selected",style={"flexBasis":"0px","flexGrow":0}),
        dcc.Store(id="resize_state",data={"width":{"is_resizing":False},"height":{"is_resizing":False}}),
    ],style={"height":"100%"},id="full_row")],style={"height":"100%"},fluid=True)
# ,style={
#         "flexGrow":1,
#         "flexShrink":0,
#         "flewBasis":"80vh",
#         #"height":"100vh",
#              "overflow":"hidden"},id="full_row")
