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
import dash  # Needed for dash.callback_context in callbacks

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
                            dbc.Tab(id="mono_tab",tab_id="mono_signature_view",label="mono signature view",children=[
                                html.Div(
                                    style={"position": "relative", "width": "100%", "height": "100%"},
                                    children=[
                                        html.Span(tooltip.create_tooltip(mono_graph),style={"position":"relative","zIndex":"4"}),
                                        mono_graph,
                                        html.Div(className="legend_canvas",children=[
                                            html.Canvas(id="mono_canvas")
                                        ])
                                    ]
                                )
                            ]),
                            dbc.Tab(id="multi_tab",tab_id="multi_signature_view",label="multi signature view",children=[
                                html.Div(
                                    style={"position": "relative", "width": "100%", "height": "100%"},
                                    children=[
                                        html.Span(tooltip.create_tooltip(multi_graph),style={"position":"relative","zIndex":"4"}),
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
                                        # Add the button in the bottom right corner, styled to fit the app, but make it conditionally visible
                                        html.Div(
                                            dbc.Button(
                                                "Show pathways",
                                                id="show_pathways_btn",
                                                color="secondary",
                                                size="sm",
                                                style={
                                                    "background": "#f8f9fa",
                                                    "color": "#222",
                                                    "border": "1px solid #ccc",
                                                    "boxShadow": "0 1px 2px rgba(0,0,0,0.04)",
                                                    "borderRadius": "6px",
                                                    "padding": "2px 10px",
                                                    "fontSize": "0.95rem",
                                                    "fontWeight": 500,
                                                    "minWidth": "0",
                                                    "minHeight": "0",
                                                }
                                            ),
                                            id="show_pathways_btn_container",
                                            style={
                                                "position": "absolute",
                                                "bottom": "-14px",
                                                "right": "5px",
                                                "zIndex": 300,
                                                "pointerEvents": "auto"
                                            }
                                        ),
                                        dcc.Store(id="legend_tooltip_store"),
                                        dcc.Store(id="legend_hover_store"),
                                        dcc.Store(id="legend_active_store"),
                                        # Add a store to track pathway visibility
                                        dcc.Store(id="show_pathways_store", data=False),
                                    ]
                                )
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
                             dcc.Store(id="pathway_signature_toast_store",data={}),
                             dcc.Store(id="pending_pathway_store", data=None),
                             dcc.Store(id="previous_selection_store", data=None),
                             dcc.Store(id="proceed_yes_store", data=None),
                             dbc.Toast(
                                 id="pathway_signature_toast",
                                 header=None,
                                 is_open=False,
                                 duration=None,  # No auto-dismiss
                                 icon=None,
                                 dismissable=False,  # Remove close (cross) button
                                 style={
                                     "position": "fixed",
                                     "top": "50%",
                                     "left": "50%",
                                     "transform": "translate(-50%, -50%)",
                                     "zIndex": 9999,
                                     "minWidth": "350px",
                                     "maxWidth": "90vw",
                                     "textAlign": "center",
                                     "background": "#fff",
                                     "border": "1px solid #e0e0e0",
                                     "boxShadow": "0 4px 24px rgba(0,0,0,0.08)",
                                     "borderRadius": "10px",
                                     "fontSize": "1.08rem",
                                     "color": "#222",
                                     "padding": "1rem 1.5rem 1rem 1.5rem"
                                 },
                                 children=[
                                     "This action will take about 1 minute, do you want to proceed?",
                                     html.Div([
                                         dbc.Button("Yes", id="proceed_yes_btn", color="primary", style={"marginRight": "1rem", "minWidth": "80px"}),
                                         dbc.Button("No", id="proceed_no_btn", color="secondary", style={"minWidth": "80px"})
                                     ], style={"marginTop": "1.2rem", "display": "flex", "justifyContent": "center", "gap": "0.5rem"})
                                 ]
                             ),
                             dbc.Toast(
                                 id="copy_toast",
                                 header="Copy to Clipboard",
                                 is_open=False,
                                 duration=3000,
                                 icon="success",
                                 style={"position": "fixed", "bottom": 20, "right": 20, "zIndex": 9999},
                                 children=""
                             ),
                             dbc.Alert(
                                 id="pathway_signature_warning_alert",
                                 is_open=False,
                                 dismissable=True,
                                 style={
                                     "position": "fixed",
                                     "top": "20px",
                                     "right": "20px",
                                     "zIndex": 9999,
                                     "maxWidth": "400px"
                                 }
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
            dcc.Store(id="revert_in_progress_store", data=False),
            dcc.Store(id="block_ui_store", data=False),
            html.Div(id="block_ui_overlay", style={"display": "none"}),
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

@callback(
    Output("pathway_signature_toast", "is_open"),
    Output("pathway_signature_toast", "children"),
    Output("pending_pathway_store", "data", allow_duplicate=True),
    Output("previous_selection_store", "data", allow_duplicate=True),
    Output("revert_in_progress_store", "data", allow_duplicate=True),
    Output("block_ui_store", "data", allow_duplicate=True),
    Input("pathway_signature_toast_store", "data"),
    Input("pathway_signature_toast", "is_open"),
    Input("proceed_no_btn", "n_clicks"),
    State("pathway_signature_toast", "is_open"),
    State("revert_in_progress_store", "data"),
    prevent_initial_call=True
)
def update_pathway_signature_toast(toast_data, toast_is_open_input, no_btn_click, toast_is_open_state, revert_in_progress):
    ctx = dash.callback_context
    triggered_revert = False
    if (
        ctx.triggered
        and (
            (ctx.triggered[0]["prop_id"] == "pathway_signature_toast.is_open" and toast_is_open_input is False and toast_is_open_state is True)
            or (ctx.triggered[0]["prop_id"] == "proceed_no_btn.n_clicks" and no_btn_click)
        )
        and not revert_in_progress
    ):
        triggered_revert = True
    if triggered_revert:
        return False, "", 'REVERT', None, True, True
    if toast_data and toast_data.get("is_open"):
        return True, toast_data.get("message", ""), dash.no_update, dash.no_update, False, False
    return dash.no_update, dash.no_update, dash.no_update, dash.no_update, False, False

@callback(
    Output("proceed_yes_store", "data", allow_duplicate=True),
    Output("pathway_signature_toast", "is_open", allow_duplicate=True),
    Input("proceed_yes_btn", "n_clicks"),
    State("pending_pathway_store", "data"),
    prevent_initial_call=True
)
def handle_proceed_yes(n_clicks, pending_pathway):
    if n_clicks and pending_pathway:
        return pending_pathway, False
    return dash.no_update, dash.no_update

@callback(
    Output("block_ui_overlay", "style"),
    Input("block_ui_store", "data"),
    prevent_initial_call=False
)
def show_block_ui_overlay(block):
    if block:
        return {
            "position": "fixed",
            "top": 0,
            "left": 0,
            "width": "100vw",
            "height": "100vh",
            "background": "rgba(255,255,255,0.6)",
            "zIndex": 99999,
            "display": "flex",
            "alignItems": "center",
            "justifyContent": "center"
        }
    else:
        return {"display": "none"}

# Add a callback at the end of the file to control the button's visibility
from dash import callback, Output, Input, State
import dash

@callback(
    Output("show_pathways_btn_container", "style"),
    Input("detail_graph", "elements"),
    State("show_pathways_btn_container", "style"),
    prevent_initial_call=False
)
def toggle_show_pathways_btn(elements, style):
    # Show the button only if the detailed graph has elements (i.e., is active and showing something)
    if elements and len(elements) > 0:
        style = dict(style) if style else {}
        style["display"] = "block"
        return style
    else:
        style = dict(style) if style else {}
        style["display"] = "none"
        return style

# Add a callback to toggle pathway visibility and update button label
from dash import callback, Output, Input, State
import dash

@callback(
    Output("show_pathways_store", "data"),
    Output("show_pathways_btn", "children"),
    Input("show_pathways_btn", "n_clicks"),
    State("show_pathways_store", "data"),
    prevent_initial_call=True
)
def toggle_pathways(n_clicks, current):
    if n_clicks is None:
        raise dash.exceptions.PreventUpdate
    new_state = not current
    return new_state, ("Hide pathways" if new_state else "Show pathways")
