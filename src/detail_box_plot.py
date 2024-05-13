from dash_extensions.enrich import dcc
import plotly.express as px
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from dash_extensions.enrich import Dash, html, Input, Output, State, callback
from collections import defaultdict

def detail_box_plot():
    fig = go.Figure(data=[])
    return dcc.Graph(figure=fig,id='activation_boxplot',config={"scrollZoom":False,"displayModeBar":False},style={"height":"96%","cursor":"auto"})

def detail_heatmap():
    fig = go.Figure(data=[])
    return dcc.Graph(figure=fig,id='activation_heatmap',config={"scrollZoom":False,"displayModeBar":False})
activation_data = pd.read_csv('data/activations/fake_data.csv',delimiter=";",dtype=defaultdict(lambda :float,Disease=str,Stage=str))
activation_data["box_category"] = activation_data["Disease"]+"_"+activation_data["Stage"]

# @callback(
#     Output("activation_boxplot","figure"),
#     Output("activation_boxplot","className"),
#     Input("overview_graph","selectedNodeData" ),
#     Input("detail_graph","selectedNodeData")
# )
# def update_box_plot(overview_selected,detail_selected):
#     if(detail_selected is not None and len(detail_selected)==1):
#         g = detail_selected[0]['id']
#         diseases = set([i['Disease'] for i in overview_selected])
#         only_selected_diseases = activation_data[activation_data["Disease"].isin(diseases)]
#         selected_patient_and_genes = only_selected_diseases.filter(["Disease","Stage",g,"box_category"])
#         box_categories = pd.unique(selected_patient_and_genes["box_category"]).tolist()
#         y = [selected_patient_and_genes[selected_patient_and_genes["box_category"]==c][g].to_numpy().tolist() for c in box_categories]
#         print("ubp",detail_selected,overview_selected,selected_patient_and_genes,box_categories,y)
#         return             px.box(selected_patient_and_genes,x=["box_category"],y=g),"visible_plot"
        
#     else:
#         return go.Figure(data=[
#         ]),"hidden_plot"

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
#         print(selected_box)
#         disease,stage = selected_box['points'][0]['x'].split("_")

#         only_selected_disease = activation_data[activation_data["Disease"]==disease]
#         selected_patient_and_genes = only_selected_disease.filter(["Disease","Stage",g,"box_category"]).sort_values(by=["Stage",g])
#     #     selected_patient_and_genes["box_category"] = selected_patient_and_genes["Disease"]+"_"+selected_patient_and_genes["Stage"]
#     #     box_categories = pd.unique(selected_patient_and_genes["box_category"]).tolist()
#     #     y = [selected_patient_and_genes[selected_patient_and_genes["box_category"]==c][g].to_numpy().tolist() for c in box_categories]
#     #     print("ubp",detail_selected,overview_selected,selected_patient_and_genes,box_categories,y)
#     #     return             px.box(selected_patient_and_genes,x=["box_category"],y=g),"visible_plot"
#     #
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

    