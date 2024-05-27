from dash_extensions.enrich import Dash, html, Input, Output, State, callback
import dash_cytoscape as cyto
from dash_extensions.enrich import clientside_callback,ClientsideFunction, Input, Output

# @callback(Output('detail_graph','elements'),
#           Input('overview_graph','selectedNodeData')
# )
# def display_detail_graph(data):
#     edges = []
#     nodes = set()
#     if data is not None:
#         for sign in data:
#             genes = sign['genes'].split(";")
#             for g in genes:
#                 nodes.add(g)
#             for i in range(len(genes)):
#                 for j in range(i+1,len(genes)):
#                     edges.append({'data':{'source':genes[i],'target':genes[j]}})
#     nodes = [{'data':{'id':g,"label":g}} for g in nodes]
#     print(data,nodes,edges)
    # return nodes+edges
