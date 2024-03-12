from dash import Dash, html, Input, Output, State, callback
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
class Controller(object):
    _instance = None
    def __new__(cls,dm):
        if cls._instance is None:
            cls._instance = super(Controller, cls).__new__(cls)
            cls._instance.dm = dm
        return cls._instance
    
@callback(Output('overview_graph','elements'),
          Input("disease_filter","value"),
          Input("stage_filter","value"))
def update_overview(diseases,stages):
    return ov.get_elements(Controller._instance.dm,disease_filter=diseases,stage_filter=stages)

def circle_approx_points(center,radius,n_points):
    angles = np.linspace(0,np.pi*2,n_points)
    return np.stack([np.cos(angles)*radius+center[0],np.sin(angles)*radius+center[1]],-1).tolist()
@callback(Output('detail_graph','elements'),
          Output('detail_graph','stylesheet'),
          Output('detail_graph','layout'),
          Input('overview_graph','selectedNodeData'),
          State('detail_graph','elements')
          )
def display_detail_graph(data,existing_elements):
        #TODO add metanodes siwth style shape + label
        def size(degree):
            default_diam = 20
            default_area = (default_diam**2)*np.pi/4
            area  = default_area*degree
            return np.sqrt(area/(np.pi/4))
            
        # updated_elements = Patch()
        edges = []
        nodes = []
        stylesheet_detail = [{
            "selector":"node",
            "style":{
                "width":"data(weight)",
                "height":"data(weight)",
                "label":"data(label)"
                }
        }]
        if data is not None and len(data)>0:
            color_by_diseases = True # TODO
            cm = Controller._instance.dm.get_disease_cmap() if color_by_diseases else Controller._instance.dm.get_stage_cmap()

            # print(updated_elements)
            # print(existing_elements)
            # updated_elements.clear()
            intersections,items = Controller._instance.dm.get_genes_intersections(id_filter=[i["id"] for i in data])

            sizes = items.set_index("gene")
            updated = set()
            print(sizes.head())
            to_del=[]
            for i in range(len(existing_elements)):
                elem = existing_elements[i]
                if not "source" in elem["data"] and elem["data"]["id"] in sizes.index:
                    elem["data"]["weight"]=size(sizes.loc[elem["data"]["id"] ]["size"])
                    elem["data"]["Signatures"]=sizes.loc[elem["data"]["id"] ]["id"]
                    elem["data"]["label"]=elem["data"]["id"] +" " +str(sizes.loc[elem["data"]["id"] ]["size"])
                    # updated_elements[i]["data"]["weight"] = size(sizes.loc[elem["data"]["id"] ]["size"])
                    # updated_elements[i]["data"]["label"] = elem["data"]["id"] +" " +str(sizes.loc[elem["data"]["id"] ]["size"])
                    updated.add(elem["data"]["id"])
                    stylesheet_detail.append({
                        "selector":f"node#{elem['data']['id']}",
                        "style":{
                            "width":f'{size(sizes.loc[elem["data"]["id"] ]["size"])}',
                            "height":f'{size(sizes.loc[elem["data"]["id"] ]["size"])}',
                            "label":"data(label)"
                            }
                    })
                else:
                    to_del.append(i)
                if elem["data"]["id"]=="AJUBA":
                    print(elem)
            nodes = [{'data':{'id':g.gene,"label":g.gene,"weight":size(g.size),"Signatures":g.id}} for g in items.itertuples() if g.gene not in updated]
            for i in reversed(to_del):
                # del updated_elements[i]
                del existing_elements[i]
            #         elem["data"]["label"]=elem["data"]["id"] + " "+str(sizes.loc[elem["data"]["id"] ]["size"])
            #         print(elem)
                
            # for g in items.itertuples():
            #     if g.gene =="AJUBA":
            #         print(g,size(g.size),)
            #     # stylesheet.append({"selector":"#"+g.gene,"style":{"width":g.size,"height":g.size}})
            edges = [{"data":{"source":k[0],"target":k[1]}} for k,v in intersections.items()]
            # updated_elements.extend(nodes)
            # updated_elements.extend(edges)

            existing_elements = existing_elements +nodes
            existing_elements = existing_elements +edges
            pathways_nodes = dict()
            pathways_stylesheets = dict()
            pathways_edges = []
            for elem in existing_elements:
                if "source" not in elem["data"]:
                    pathways = Controller._instance.dm.get_pathways([elem['data']['id']])
                    for p in pathways[elem['data']['id']]:
                        if not p in pathways_nodes:
                            pathways_nodes[p]={"data":{"id":p,"label":p,"weight":1},"selectable":False}
                            pathways_stylesheets[p]={
                                "selector":f"node#{p}",
                                "style":{
                                    "shape":"hexagon",
                                    # "background-opacity":0.25,
                                    # "font-size":5,
                                    "z-index":-1,
                                    # "text-halign":"left",
                                    # "text-valign":"center",
                                    # "border-width":1,
                                    # "text-wrap":"wrap",
                                    "label":p
                                    }
                                }
                        else:
                            pathways_nodes[p]["data"]["weight"]=pathways_nodes[p]["data"]["weight"]+1
                        pathways_edges.append({"data":{"source":elem["data"]["id"],"target":p}})
            for p in pathways_nodes:
                pathways_nodes[p]["data"]["weight"]=size(pathways_nodes[p]["data"]["weight"])

            existing_elements = existing_elements + list(pathways_nodes.values()) + pathways_edges
            stylesheet_detail  = stylesheet_detail + list(pathways_stylesheets.values())

            pos = nx_layout.spring_layout(existing_elements)
            points = {}

            for elem in existing_elements:
                if "source" not in elem["data"]:
                    p = pos[elem["data"]['id']]
                    if "position" in elem:
                        elem["position"]["x"]=p[0]
                        elem["position"]["y"]=p[1]
                    else:
                        elem["position"]={'x':p[0],'y':p[1]}
                    # added_points =     [[p[0]-elem["data"]['weight'],p[1]-elem["data"]['weight']],
                                # [p[0]+elem["data"]['weight'],p[1]-elem["data"]['weight']],
                                # [p[0]-elem["data"]['weight'],p[1]+elem["data"]['weight']],
                                # [p[0]+elem["data"]['weight'],p[1]+elem["data"]['weight']]]
                    added_points = circle_approx_points(p,elem["data"]['weight'],12)
                    if  "Signatures" in elem["data"]:
                        for i in elem["data"]["Signatures"]:
                            if i in points:
                                points[i] = points[i] + added_points
                            else:
                                points[i] = added_points

                            
            for i in points:
                sign_pos = np.array(points[i])
                hull = ConvexHull(sign_pos)
                bb = (np.min(sign_pos[hull.vertices],axis=0),np.max(sign_pos[hull.vertices],axis=0))
                shape_polygon_points  = []

                for v in hull.vertices:
                    shape_polygon_points.append(-1 + 2*(sign_pos[v,0]-bb[0][0])/(bb[1][0]-bb[0][0]))
                    shape_polygon_points.append(-1 + 2*(sign_pos[v,1]-bb[0][1])/(bb[1][1]-bb[0][1]))
                # print(i,sign_pos,hull.vertices,shape_polygon_points)
                existing_elements.append({
                    "data":{"id":i,"label":i},
                    "position":{"x":0.5*(bb[0][0]+bb[1][0]),
                                "y":0.5*(bb[0][1]+bb[1][1])},
                    "selectable": False,
                    "classes":" ".join([j for j in i.split("_") if j!="vs"])})
                stylesheet_detail.append({
                    "selector":f"node#{i}",
                    "style":{
                        "width":f"{bb[1][0]-bb[0][0]}",
                        "height":f"{bb[1][1]-bb[0][1]}",
                        "shape":"polygon",
                        "background-opacity":0.25,
                        "font-size":5,
                        "z-index":-1,
                        "text-halign":"left",
                        "text-valign":"center",
                        "border-width":1,
                        "text-wrap":"wrap",
                        "shape-polygon-points":" ".join(list(map(str,shape_polygon_points))),
                        "label":"\n".join(i.split("_"))
                        }
                     })
            for d in cm:
                color_str = f"rgba({','.join([str(i*255) for i in cm[d][:3]])},0.25)"

                stylesheet_detail.append({
                    "selector":f"node.{d}",
                    "style":{
                        "background-color":color_str
                    }
                })
            # shape polygon
                        # shape-polygon-points
                        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
            # updated_elements.extend(edges)

        # if data is not None:
            # nodes = set()
        #     for sign in data:
        #         genes = sign['genes'].split(";")
        #         for g in genes:
        #             nodes.add(g)
        #         for i in range(len(genes)):
        #             for j in range(i+1,len(genes)):
        #                 edges.append({'data':{'source':genes[i],'target':genes[j]}})
        # nodes = [{'data':{'id':g,"label":g}} for g in nodes]
        # print(nodes)
        return existing_elements,stylesheet_detail,{"name":"preset",'positions': {
                node['data']['id']: node['position']
                for node in existing_elements if "source" not in node["data"]
            }}



@callback(
    Output("activation_boxplot","figure"),
    Output("activation_boxplot","className"),
    Output("overview_graph","stylesheet"),
    State("overview_graph","selectedNodeData" ),
    Input("detail_graph","selectedNodeData"),
        State("overview_graph","elements"),
        State("overview_graph","stylesheet"),

)
def update_box_plot(overview_selected,detail_selected,overview_elements,overview_stylesheets):
    if(detail_selected is not None and len(detail_selected)>=1):
        g = detail_selected[0]['id']
        items = [i["id"] for i in detail_selected]
        diseases = set([i['Disease'] for i in overview_selected])
        selected_patient_and_genes = Controller._instance.dm.get_activations(g,diseases)
        box_categories = pd.unique(selected_patient_and_genes["box_category"]).tolist()
        # print([ i for i in overview_elements if "elems" in i["data"] and  any([j in items for j in i["data"]["elems"]])])
        stylesheets =overview_stylesheets if overview_stylesheets is not None else ov.get_default_stylesheet(Controller._instance.dm)
        stylesheets = [s for s in stylesheets if not(s["selector"].startswith("edge#"))]
        for i in overview_elements:
            if "elems" in i["data"] and  any([j in items for j in i["data"]["elems"]]):

                stylesheets.append({"selector":'edge#'+i["data"]["id"],"style":{"line-color":"red"}})
                # if("classes" in i):
                    # i["classes"] = " ".join(set(i["classes"].split(" ") + ["highlight_edge"]))
                # else:
                    # i["classes"] =  ["highlight_edge"]
            # else:
                # if("classes" in i) and "highlight_edge" in i["classes"]:
                    # classes:list = i["classes"].split(" ")
                    # classes.remove("highlight_edge")
                    # i["classes"] = " ".join(classes)
        return px.box(selected_patient_and_genes,x=["box_category"],y=g),"visible_plot",stylesheets
        
    else:
        return go.Figure(data=[
            ]),"hidden_plot",overview_stylesheets if overview_stylesheets is not None else ov.get_default_stylesheet(Controller._instance.dm)
    
@callback(
    Output("activation_heatmap","figure"),
    Output("activation_heatmap","className"),
    Input("overview_graph","selectedNodeData" ),
    Input("detail_graph","selectedNodeData"),
    Input("activation_boxplot","clickData")
)
def update_heatmap(overview_selected,detail_selected,selected_box):
    if(detail_selected is not None and len(detail_selected)==1 and selected_box is not None and len(selected_box)==1):
        g = detail_selected[0]['id']
        # print(selected_box)
        disease,stage = selected_box['points'][0]['x'].split("_")

        selected_patient_and_genes = Controller._instance.dm.get_activations(g,[disease]).sort_values(by=["Stage",g])
        y_labels = []
        for item in selected_patient_and_genes.itertuples():
            if len(y_labels)==0 or y_labels[-1]!=item.box_category:
                y_labels.append(item.box_category)
            else:
                y_labels.append("")
        fig =  px.imshow(np.expand_dims(selected_patient_and_genes[g].to_numpy(),1),y=y_labels,color_continuous_scale="turbo")
        fig.update_xaxes(showticklabels=False)

        return fig,"visible_plot"
    else:
        return go.Figure(data=[
        ]),"hidden_plot"

    