from dash import Dash, html, Input, Output, State, callback,ctx
import dash_cytoscape as cyto
from dash import clientside_callback,ClientsideFunction, Input, Output,Patch
import numpy as np
import plotly.graph_objs as go
import pandas as pd
import plotly.express as px
import overview as ov
import nx_layout
import scipy
from scipy.spatial import ConvexHull
import detail_graph
import utils
class Controller(object):
    _instance = None
    def __new__(cls,dm):
        if cls._instance is None:
            cls._instance = super(Controller, cls).__new__(cls)
            cls._instance.dm = dm
        return cls._instance

@callback(
        Output('overview_graph','elements'),
          Input("disease_filter","value"),
          Input("comparisons_filter","value"),
          Input("overview_graph","selectedNodeData"),
          Input("detail_graph","selectedNodeData"),
          Input("genes_menu_select","value"),
          State("overview_graph","elements")
          )
def update_overview(diseases,comparisons_filter,selectedSign,selected_detail_gene,menu_selected_genes,cur_elems):
    
    classes = ["highlight","half_highlight"]
    for i in cur_elems:
        if "Disease" in i["data"]:
            if i['data']["Disease"] in diseases and i['data']["Comparison"] in comparisons_filter:
                c = "highlight"
                if( menu_selected_genes is not None and menu_selected_genes != "None"):
                    c = "half_highlight"
                    if menu_selected_genes in i["data"]["Signature"]:
                        c = "highlight"
                i["classes"]=utils.switch_class(i["classes"],[c],classes)
            else:
                i["classes"]=utils.switch_class(i["classes"],[],classes)
        else:
            src = i["data"]["source"].split("_")
            tgt = i["data"]["target"].split("_")
            src.remove("vs")
            tgt.remove("vs")
            src_dis = src[0]
            src_comp = "_".join(src[1:])
            tgt_comp = "_".join(tgt[1:])
            tgt_dis = tgt[0]
            if(src_dis in diseases and tgt_dis in diseases and src_comp in comparisons_filter and tgt_comp in comparisons_filter):
                if "highlight" not in i["classes"]:
                    i["classes"] = " ".join(i["classes"].split(" ")+["highlight"])
            else:
                if "highlight" in i["classes"]:
                    classes = i["classes"].split(" ")
                    classes.remove("highlight")
                    i["classes"]="_".join(classes)
    return cur_elems

@callback(Output('data_overview_selected','value'),
          Input('overview_graph','selectedNodeData'), prevent_initial_call=True
          )
def update_overview_selected_data(data):
    if data is None:
        return ""
    return ";".join([i["id"] for i in data])
@callback(Output('data_menu_selected','value'),
          Input('disease_filter','value'),
            prevent_initial_call=True
          )
def update_menu_selected_data(data):
    if data is None:
        return ""
    return ";".join(data)
@callback(Output("selected_genes","data"),
          Input("genes_menu_select","value"))
def select_genes(menu_select):
    if menu_select is not None and menu_select != "None":
        return [menu_select]
    else:
        return []
@callback(Output('disease_filter','value'),
          Input("selected_genes","data"),
            prevent_initial_call=True
          )
def update_menu_disease_filter_data(genes):
    return Controller._instance.dm.get_diseases_from_genes(genes)
@callback(
        Output('comparisons_filter','value'),
          Input("selected_genes","data"),
        State('comparisons_filter','value'),
            prevent_initial_call=True
          )
def update_menu_comparison_filter_data(genes,cur_state):
    if len(genes)>0:
        return Controller._instance.dm.get_comparisons_from_genes(genes)
    return cur_state

@callback(Output('detail_graph','elements'),
          Output('detail_graph','stylesheet'),
          Output('detail_graph','layout'),
        #   Input('overview_graph','selectedNodeData'),
          Input('disease_filter','value'),
          Input('comparisons_filter','value'),
          Input('data_overview_selected','value'),
          Input('selected_genes','data'),
          State('detail_graph','elements'),
            prevent_initial_call=True
          )
def display_detail_graph(diseases,comparisons,signatures,menu_genes,existing_elements):
    if(all([len(diseases)==0 or len(comparisons)==0 ,signatures is None])):
        return [],[],{"name":"preset"}
    if len(diseases)!=0 or len(signatures)!=0:
        if diseases is None:
            diseases =""
        if signatures is None:
            signatures = ""
        return detail_graph.display_detail_graph(list(filter(lambda a: len(a)>0,diseases)),list(filter(lambda a: len(a)>0,signatures.split(";"))),menu_genes,existing_elements)
    else:
        return existing_elements,[],{"name":"preset"}

@callback(Output("data_gene_detail_selected","value"),
          Input("detail_graph","selectedNodeData"))
def data_gene_detail_selected(nodes):
    if nodes is None:
        return ""
    return ";".join([i["id"] for i in nodes])

@callback(
    Output("activation_boxplot","figure"),
    Output("activation_boxplot","className"),
    Output("overview_graph","stylesheet"),
    Input("data_menu_selected","value"),
    Input("data_overview_selected","value"),
    Input("data_gene_detail_selected","value"),
    # State("overview_graph","selectedNodeData" ),
    Input("selected_genes","data"),
        State("overview_graph","elements"),
        State("overview_graph","stylesheet"),
        prevent_initial_call=True
)
def update_box_plot(menu_selected_diseases,overview_selected,detail_selected,menu_selected,overview_elements,overview_stylesheets):
    items = []
    if(detail_selected is not None and len(detail_selected)>=1):
        items = detail_selected.split(";")
    elif menu_selected is not None and len(menu_selected)!=0:
        items = menu_selected
    if(len(items)>0):
        g = items[0]
        diseases = menu_selected_diseases.split(";")
        if overview_selected is not None:
            diseases = diseases + list(set([i['Disease'] for i in overview_selected]))
        selected_patient_and_genes = Controller._instance.dm.get_activations(items,diseases)
        box_categories = pd.unique(selected_patient_and_genes["box_category"]).tolist()
        stylesheets =overview_stylesheets if overview_stylesheets is not None else ov.get_default_stylesheet(Controller._instance.dm)
        stylesheets = [s for s in stylesheets if not(s["selector"].startswith("edge#"))]
        for i in overview_elements:
            if "elems" in i["data"] and  any([j in items for j in i["data"]["elems"]]):
                stylesheets.append({"selector":'edge#'+i["data"]["id"],"style":{"line-color":"red","width":6}})
                # if("classes" in i):
                    # i["classes"] = " ".join(set(i["classes"].split(" ") + ["highlight_edge"]))
                # else:
                    # i["classes"] =  ["highlight_edge"]
            # else:
                # if("classes" in i) and "highlight_edge" in i["classes"]:
                    # classes:list = i["classes"].split(" ")
                    # classes.remove("highlight_edge")
                    # i["classes"] = " ".join(classes)
        if(len(items)==1):
            return px.box(selected_patient_and_genes,x="box_category",y=items[0]),"visible_plot",stylesheets
        else:
            dfs = []
            for i in range(len(items)):
                df = selected_patient_and_genes.filter(("box_category",items[i]))
                df = df.rename({items[i]:"expression"},axis=1)
                df["gene"] = items[i]
                dfs.append(df)
            df = pd.concat(dfs)
            return px.box(df,x="box_category",y="expression",color="gene"),"visible_plot",stylesheets
        
    else:
        return go.Figure(data=[
            ]),"hidden_plot",overview_stylesheets if overview_stylesheets is not None else ov.get_default_stylesheet(Controller._instance.dm)
    
# @callback(
#     Output("activation_heatmap","figure"),
#     Output("activation_heatmap","className"),
#     Input("overview_graph","selectedNodeData" ),
#     Input("detail_graph","selectedNodeData"),
#     Input("activation_boxplot","clickData")
# )
# def update_heatmap(overview_selected,detail_selected,selected_box):
#     if(detail_selected is not None and len(detail_selected)==1 and selected_box is not None and len(selected_box)==1):
#         g = detail_selected[0]['id']
#         disease,stage = selected_box['points'][0]['x'].split("_")

#         selected_patient_and_genes = Controller._instance.dm.get_activations(g,[disease]).sort_values(by=["Stage",g])
#         y_labels = []
#         for item in selected_patient_and_genes.itertuples():
#             if len(y_labels)==0 or y_labels[-1]!=item.box_category:
#                 y_labels.append(item.box_category)
#             else:
#                 y_labels.append("")
#         fig =  px.imshow(np.expand_dims(selected_patient_and_genes[g].to_numpy(),1),y=y_labels,color_continuous_scale="turbo")
#         fig.update_xaxes(showticklabels=False)

#         return fig,"visible_plot"
#     else:
#         return go.Figure(data=[
#         ]),"hidden_plot"

@callback(
        Output("data_gene_menu_selected","value"),
        Input("genes_menu_select","value")
)
def update_data_gene_menu_selected(v):
    return v


clientside_callback(
        """
function update_tooltip_data(mouseOver){
if(mouseOver!==undefined){
console.log("tooltip","not unedfined",mouseOver);
var elem = document.getElementById('detail_graph');

    return {'data':mouseOver,'width':elem.clientWidth,'height':elem.clientHeight};
}
console.log("tooltip","unedfined",mouseOver);
return {};
}
        """,
            Output("fake_graph_size","data"),
            Input("detail_graph","mouseoverNodeData"),
            )

@callback(
    Output("detail_graph_tooltip","children"),
    Output("detail_graph_tooltip","show"),
    Output("detail_graph_tooltip","bbox"),
    Output("detail_graph_tooltip","direction"),
#     # Output("detail_graph_tooltip","bbox"),
    Input("fake_graph_size","data"),
    State("detail_graph","elements"),
    State("detail_graph","extent"),
    State("detail_graph","stylesheet"),prevent_initial_call=True
    )
def activate_tooltip(mouseoverNodeData,elements,extent,stylesheets):
    if mouseoverNodeData is None or "data" not in mouseoverNodeData or mouseoverNodeData["data"] is None:
        return html.Span(""),False,{},"right"
    else:
        elem = None
        for i in elements:
            if i["data"]["id"] == mouseoverNodeData["data"]["id"]:
                elem = i
        width = None
        height =None
        wratio = (mouseoverNodeData["width"])/ (extent["x2"]-extent["x1"])
        hratio = (mouseoverNodeData["height"])/ (extent["y2"]-extent["y1"])

        x = (elem["position"]["x"]-extent["x1"]) *wratio
        y = (elem["position"]["y"]-extent["y1"]) *hratio
        if(x/mouseoverNodeData["width"]>0.5):
            direction="left"
        else:
            direction="right"
        if "weight" in mouseoverNodeData["data"]:
            width = mouseoverNodeData["data"]["weight"]
            height = mouseoverNodeData["data"]["weight"]
        else:
            for i in stylesheets:
                if i["selector"] == f'node#{mouseoverNodeData["data"]["id"]}':
                    width = float(i["style"]["width"])
                    height = float(i["style"]["height"])
                    polygon = np.array(list(map(float,i["style"]["shape-polygon-points"].split(" ")))).reshape(-1,2)
            zeros = np.argwhere(polygon[:,1]==0)
            y_crosses = []
            for z in zeros:
                y_crosses.append(polygon[z,0])
            if (len(y_crosses)<2):
                for i in range(len(polygon)):
                    if(polygon[i,1]!=0.0 and polygon[(i+1)%len(polygon),1]!=0.0):
                        if(polygon[i,1] * polygon[(i+1)%len(polygon),1]<0):
                            y_crosses.append(
                                polygon[i,0]+(abs(polygon[i,1])/(abs(polygon[i,1])+abs(polygon[(i+1)%len(polygon),1])))*(polygon[(i+1)%len(polygon),0]-polygon[i,0])
                                )
            if(len(y_crosses)==2):
                if(direction=="right"):
                    x+=0.5*np.max(y_crosses)*width*wratio
                else:
                    x+=0.5*np.min(y_crosses)*width*wratio
                width=0
                height=0
                


        return html.Span(mouseoverNodeData["data"]["id"]),True,{"x0":x-width*wratio/2,"y0":y-height*hratio/2,"x1":x+width*wratio/2,"y1":y+height*hratio/2},direction
    
# clientside_callback(ClientsideFunction(
#     namespace="tooltip",
#     function_name="set_tooltip"),
#     Output("dummy_div","children"),
#     Input("detail_graph","mouseoverNodeData"),prevent_initial_call=True
# )