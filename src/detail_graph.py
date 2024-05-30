import re
import time
import dash_cytoscape as cyto
from dash import Patch
import numpy as np
import pandas as pd
import nx_layout
import tlp_layout
from controller import Controller
from scipy.spatial import ConvexHull
import plotly.express as px
from shapely import Polygon,Point
from dash_extensions.enrich import html
from dash_svg import Svg, G, Path, Circle,Rect
# import euler_layout
class Edge:
    def __init__(self,source_id,target_id):
        self.source_id = source_id
        self.target_id = target_id
    def __eq__(self,other):
        return self.source_id==other.source_id and self.target_id == other.target_id
    def __hash__(self):
        return (hash(self.source_id)<<32)^hash(self.target_id)
    @property
    def data(self):
        return {'data':{'source':self.source_id,'target':self.target_id}}
# app.clientside_callback(
# """function (i) {
#     if (i) {
#         yourfunction()
# }
# return window.dash_clientside.no_update
# }""",
# Ouput('buttonid','id'),
# Input('buttonid','n_clicks')
# )
def detail_graph():
    return cyto.Cytoscape(
        id="detail_graph",
        style={
            "width":"97%",
            "height":"100%"},
        autoungrabify=True,
                wheelSensitivity=0.5,
                maxZoom=4,
                minZoom=0.1,
        layout={"name":"preset"},
        clearOnUnhover=True,
        className="detail_graph",stylesheet=[{"selector":"node","style":{"label":"nothing"}}],
        elements=
            []
    )


def size(degree):
    default_diam = 20
    default_area = (default_diam**2)*np.pi/4
    area  = default_area*degree
    return np.sqrt(area/(np.pi/4))
def size_triangle(degree):
    default_height = 20
    default_area = (default_height**2)/2
    area  = default_area*degree
    return np.sqrt(2*area)

def circle_approx_points(center,radius,n_points):
    angles = np.linspace(0,np.pi*2,n_points)
    return np.stack([np.cos(angles)*radius*(1+np.random.uniform(0.05,0.15))+center[0],np.sin(angles)*radius*(1+np.random.uniform(0.05,0.15))+center[1]],-1).tolist()

def update_cur_elements(existing_elements,updated,stylesheet_detail,sizes,updated_elements):
    to_del=[]
    for i in range(len(existing_elements)):
        elem = existing_elements[i]
        if not "source" in elem["data"] and elem["data"]["id"] in sizes.index:
            elem["data"]["weight"]=size(sizes.loc[elem["data"]["id"] ]["size"])
            elem["data"]["Signatures"]=sizes.loc[elem["data"]["id"] ]["id"]
            updated_elements[i]["data"]["weight"]=size(sizes.loc[elem["data"]["id"] ]["size"])
            updated_elements[i]["data"]["Signatures"]=sizes.loc[elem["data"]["id"] ]["id"]
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
    for i in reversed(to_del):
        del existing_elements[i]
        del updated_elements[i]    

def add_path_ways(existing_elements,stylesheet_detail,updated_elements):
    pathways_nodes = dict()
    protein_nodes = dict()
    pathways_stylesheets = dict()
    pathways_edges = []
    pathways = Controller._instance.dm.get_pathways([elem['data']['id'] for elem in existing_elements if "source" not in elem["data"]])
    for p in pathways.itertuples():
        if not p.PathwayStId in pathways_nodes:
            children = []
            children.append(html.H6(p.PathwayDisplayName))
            children.append(html.A("PathwayReactomeLink",href=p.PathwayReactomeLink,target="_blank"))
            pathways_nodes[p.PathwayStId]={"data":{"id":p.PathwayStId,"label":"","weight":1,"tooltip_content":children,"is_pathways":True},"selectable":False}
            pathways_stylesheets[p.PathwayStId] = {
                        "selector":f"node#{p.PathwayStId}",
                        "style":{
                            "shape":"triangle",
                            # "background-opacity":0.25,
                            # "font-size":5,
                            "backgroundColor":"gray",

                            "z-index":-1,
                            # "text-halign":"left",
                            # "text-valign":"center",
                            # "border-width":1,
                            # "text-wrap":"wrap",
                            # "label":p.PathwayStId,
                            "width":size_triangle(1),
                            "height":size_triangle(1)
                            }
                        }
        else:
            pathways_nodes[p.PathwayStId]["data"]["weight"]=pathways_nodes[p.PathwayStId]["data"]["weight"]+1
            pathways_stylesheets[p.PathwayStId]["style"]["width"]=size_triangle(pathways_nodes[p.PathwayStId]["data"]["weight"])
            pathways_stylesheets[p.PathwayStId]["style"]["height"]=size_triangle(pathways_nodes[p.PathwayStId]["data"]["weight"])
        UniProtID = p.UniProtID 
        if UniProtID == "" or UniProtID is None or UniProtID=="None":
            UniProtID = f"prot_{p.EnsemblID}"
        if not UniProtID in protein_nodes:
            children = []
            children.append(html.H6(UniProtID))
            protein_nodes[UniProtID]={"data":{"id":UniProtID,"label":"","weight":1,"tooltip_content":children,"is_pathways":True,"t1":str(p.UniProtID),"t2":str(type(UniProtID)),"t3":str(p)},"selectable":False}
            pathways_stylesheets[UniProtID] = {
                        "selector":f"node#{UniProtID}",
                        "style":{
                            "shape":"barrel",
                            # "background-opacity":0.25,
                            # "font-size":5,
                            "backgroundColor":"gray",

                            "z-index":-1,
                            # "text-halign":"left",
                            # "text-valign":"center",
                            # "border-width":1,
                            # "text-wrap":"wrap",
                            # "label":p.PathwayStId,
                            "width":size_triangle(1),
                            "height":size_triangle(1)
                            }
                        }
        else:
            protein_nodes[UniProtID]["data"]["weight"]=protein_nodes[UniProtID]["data"]["weight"]+1
            pathways_stylesheets[UniProtID]["style"]["width"]=size_triangle(protein_nodes[UniProtID]["data"]["weight"])
            pathways_stylesheets[UniProtID]["style"]["height"]=size_triangle(protein_nodes[UniProtID]["data"]["weight"])
        
        pathways_edges.append({"data":{"source":p.EnsemblID,"target":UniProtID,"is_pathway_edge":True}})
        pathways_edges.append({"data":{"source":UniProtID,"target":p.PathwayStId,"is_pathway_edge":True}})

    # for i,elem in enumerate(existing_elements):
    #     if "source" not in elem["data"]:
    #         for p in pathways[elem['data']['id']]:
    #             if not p["kegg_id"] in pathways_nodes:
    #                 children = []
    #                 if "name" in p:
    #                     children.append(html.H6(p["name"]))
    #                 else:
    #                     print(p.keys(),p["kegg_id"])
    #                 if "description" in p:
    #                     children.append(html.P(p["description"]))
    #                 if("ko_pathway" in p):
    #                     children.append(html.A("KEGG",href=f"https://www.kegg.jp/entry/{p['ko_pathway']}",target="_blank"))
    #                 pathways_nodes[p["kegg_id"]]={"data":{"id":p["kegg_id"],"label":p["kegg_id"],"weight":1,"tooltip_content":children,"is_pathways":True},"selectable":False}
    #                 pathways_stylesheets[p["kegg_id"]] = {
    #                     "selector":f"node#{p['kegg_id']}",
    #                     "style":{
    #                         "shape":"triangle",
    #                         # "background-opacity":0.25,
    #                         # "font-size":5,
    #                         "backgroundColor":"gray",

    #                         "z-index":-1,
    #                         # "text-halign":"left",
    #                         # "text-valign":"center",
    #                         # "border-width":1,
    #                         # "text-wrap":"wrap",
    #                         "label":p["kegg_id"],
    #                         "width":size(1),
    #                         "height":size(1)
    #                         }
    #                     }
    #             else:
    #                 pathways_nodes[p["kegg_id"]]["data"]["weight"]=pathways_nodes[p["kegg_id"]]["data"]["weight"]+1
    #                 pathways_stylesheets[p["kegg_id"]]["style"]["width"]=size(pathways_nodes[p["kegg_id"]]["data"]["weight"])
    #                 pathways_stylesheets[p["kegg_id"]]["style"]["height"]=size(pathways_nodes[p["kegg_id"]]["data"]["weight"])
    #             pathways_edges.append({"data":{"source":elem["data"]["id"],"target":p["kegg_id"],"is_pathway_edge":True}})
    # pathways_edges = [e for e in pathways_edges if pathways_nodes[e["data"]["target"]]["data"]["weight"]>1]
    # pathways_nodes = {k:v for k,v in pathways_nodes.items() if v["data"]["weight"]>1}
    for p in pathways_nodes:
        pathways_nodes[p]["data"]["weight"]=size_triangle(pathways_nodes[p]["data"]["weight"])
    for elem in list(pathways_nodes.values())+ list(protein_nodes.values()) + pathways_edges:
        updated_elements.append(elem)
    return existing_elements + list(pathways_nodes.values())  + list(protein_nodes.values())+ pathways_edges,stylesheet_detail + list(pathways_stylesheets.values())


def get_genes_hull_points(existing_elements,pos,updated_elements):
    genes_hull_points = {}
    genes_anchor_points = {}

    for i,elem in enumerate(existing_elements):
        if "source" not in elem["data"] and ("is_metanode" not in elem["data"] or not elem["data"]["is_metanode"]):
            p = pos[elem["data"]['id']]
            if "position" in elem:
                elem["position"]["x"]=p[0]
                elem["position"]["y"]=p[1]
                updated_elements[i]["position"]["x"]=p[0]
                updated_elements[i]["position"]["y"]=p[1]
            else:
                elem["position"]={'x':p[0],'y':p[1]}
                updated_elements[i]["position"]={'x':p[0],'y':p[1]}
            elem["data"]["pos"]=elem["position"]
            updated_elements[i]["data"]["pos"]=elem["position"]
            if  "Signatures" in elem["data"]:
                for i in elem["data"]["Signatures"]:
                    added_points = circle_approx_points(p,elem["data"]['weight'],12)
                    if i in genes_hull_points:
                        genes_hull_points[i] = genes_hull_points[i] + added_points
                        genes_anchor_points[i] = genes_anchor_points[i] + [p]
                    else:
                        genes_hull_points[i] = added_points
                        genes_anchor_points[i] = [p]
    return genes_hull_points,genes_anchor_points

def choose_signature_label_pos(shape_polygon_points):
    hull = np.array(shape_polygon_points).reshape(-1,2)
    p = Polygon(hull)
    d = []
    for x in [-1,0,1]:
        for y in [-1,0,1]:
            if x!=0 or y!=0:
                d.append({"x":x,"y":y,"d":Point((x,y)).distance(p)})
    d_min = min(d,key=lambda x:x["d"])
    x = "left"
    y="center"
    match d_min["x"]:
        case -1:
            x = "left"
        case 0:
            x = "center"
        case 1:
            x = "right"
    match d_min["y"]:
        case 1:
            y = "bottom"
        case 0:
            y = "center"
        case -1:
            y = "top"
    return x,y
def add_signature_metanodes(gene_hull_points,existing_elements,stylesheet_detail,updated_elements,genes_anchor_points,style_only=False):
    
    index = dict([(elem["selector"],i) for i,elem in enumerate(stylesheet_detail)])

    for i in gene_hull_points:
        sign_pos = np.array(gene_hull_points[i])
        anchor_pos = np.array(genes_anchor_points[i])
        hull = ConvexHull(sign_pos)
        bb = (np.min(sign_pos[hull.vertices],axis=0),np.max(sign_pos[hull.vertices],axis=0))
        shape_polygon_points  = []
        shape_polygon_anchors  = []

        for v in hull.vertices:
            shape_polygon_points.append(-1 + 2*(sign_pos[v,0]-bb[0][0])/(bb[1][0]-bb[0][0]))
            shape_polygon_points.append(-1 + 2*(sign_pos[v,1]-bb[0][1])/(bb[1][1]-bb[0][1]))
            shape_polygon_anchors.append(-1 + 2*(anchor_pos[v//(len(sign_pos)//len(anchor_pos)),0]-bb[0][0])/(bb[1][0]-bb[0][0]))
            shape_polygon_anchors.append(-1 + 2*(anchor_pos[v//(len(sign_pos)//len(anchor_pos)),1]-bb[0][1])/(bb[1][1]-bb[0][1]))
        label_pos = choose_signature_label_pos(shape_polygon_points)
        stylesheet = {
            "selector":f"node#{i}",
            "style":{
                "width":f"{bb[1][0]-bb[0][0]}",
                "height":f"{bb[1][1]-bb[0][1]}",
                "shape":"polygon",
                "backgroundOpacity":0.25,
                "fontSize":10,
                "zIndex":1,
                "textHalign":label_pos[0],
                "textValign":label_pos[1],
                "borderWidth":1,
                "textWrap":"wrap",
                "shapePolygonPoints":" ".join(list(map(str,shape_polygon_points))),
                "label":"\n".join(i.split("_"))
                }
                }
        if not style_only:    
            existing_elements.append({
                "data":{"id":i,"label":i,"is_metanode":True,"tooltip_content":i},
                "position":{"x":0.5*(bb[0][0]+bb[1][0]),
                            "y":0.5*(bb[0][1]+bb[1][1])},
                "selectable": False,
                "classes":" ".join([j for j in i.split("_") if j!="vs"])})
            updated_elements.append({
                "data":{"id":i,"label":i,"is_metanode":True,"tooltip_content":i},
                "position":{"x":0.5*(bb[0][0]+bb[1][0]),
                            "y":0.5*(bb[0][1]+bb[1][1])},
                "selectable": False,
                "classes":" ".join([j for j in i.split("_") if j!="vs"])})
            stylesheet_detail.append(stylesheet)
        else:
            if not stylesheet["selector"] in index:
                stylesheet_detail.append(stylesheet)
            else:
                stylesheet_detail[index[stylesheet["selector"]]]=stylesheet
            updated_elements[i]["position"]={"x":0.5*(bb[0][0]+bb[1][0]),"y":0.5*(bb[0][1]+bb[1][1])}
        

def color_metanodes(cm,stylesheet_detail):
    for d in cm:
        color_str = f"rgba({','.join([str(i*255) for i in cm[d][:3]])},0.25)"

        stylesheet_detail.append({
            "selector":f"node.{d}",
            "style":{
                "background-color":color_str
            }
        })

def get_color_scheme(selected_genes):
    if len(selected_genes)<len(px.colors.qualitative.D3):
        return px.colors.qualitative.D3
    if len(selected_genes)<len(px.colors.qualitative.Safe):
        return px.colors.qualitative.Safe
    if len(selected_genes)<len(px.colors.qualitative.Dark24):
        return px.colors.qualitative.Dark24
    return px.colors.qualitative.Alphabet
def color_selected_node(stylesheet_detail,selected_genes):
    color_scheme = get_color_scheme(selected_genes)
    for i,g in enumerate(selected_genes):
        stylesheet_detail.append({
            "selector":"node#"+g,
            "style":{
                "background-color":color_scheme[i%len(color_scheme)]
                }
        })
    return stylesheet_detail

def display_detail_graph(selectedDiseases,selected_signatures,selected_genes,existing_elements,detail_pos_store,AR):
    updated_elements = Patch()
    edges = []
    nodes = []
    stylesheet_detail = [{
        "selector":"node",
        "style":{
            "width":"data(weight)",
            "height":"data(weight)",
            "label":"data(label)",
            "backgroundColor":"black",
            "zIndex":2
            }
    }]
    stylesheet_detail = color_selected_node(stylesheet_detail,selected_genes)
    color_by_diseases = True # TODO
    cm = Controller._instance.dm.get_disease_cmap()# if color_by_diseases else Controller._instance.dm.get_stage_cmap()

    # print(updated_elements)
    # print(existing_elements)
    # updated_elements.clear()

    intersections,items = Controller._instance.dm.get_genes_intersections(id_filter=selected_signatures,disease_filter=selectedDiseases,gene_filter=selected_genes)
    data = Controller._instance.dm.get_genes_intersection_data(id_filter=selected_signatures,disease_filter=selectedDiseases,gene_filter=selected_genes)
    # euler_layout.euler_layout(data)

    color_metanodes(cm,stylesheet_detail)
    if(items is None):
        return redraw(existing_elements,detail_pos_store,AR,stylesheet_detail)
    sizes = items.set_index("gene")
    updated = set()
    update_cur_elements(existing_elements,updated,stylesheet_detail,sizes,updated_elements)

        
    nodes = [{'data':{'id':g.gene,"label":"","weight":size(g.size),"Signatures":g.id,"tooltip_content":Controller._instance.dm.get_symbol(g.gene)}} for g in items.itertuples() if g.gene not in updated]

    edges = [{"data":{"source":k[0],"target":k[1]}} for k,v in intersections.items() ]

    # print(edges)
    existing_elements = existing_elements +nodes
    for n in nodes:
        updated_elements.append(n)
    # existing_elements= existing_elements +edges
    # for e in edges:
    #    updated_elements.append(e)


    existing_elements,stylesheet_detail = add_path_ways(existing_elements,stylesheet_detail,updated_elements)
    tmp_nodes = dict()
    tmp_edges = []
    for n in existing_elements+nodes:
        if "source" not in n["data"] and not ("is_pathways" in n["data"] and n["data"]["is_pathways"]):
            if "Signatures" not in n["data"]:
                print(n)
            for i in n["data"]["Signatures"]:
                if i not in tmp_nodes:
                    tmp_nodes[i]={'data':{'id':i,"weight":1}}
                tmp_edges.append({'data':{'source':n['data']['id'],"target":i,"weight":1}})
    
    pos = tlp_layout.sgd2_layout(existing_elements+list(tmp_nodes.values())+tmp_edges,detail_pos_store,AR)
    genes_hull_points,genes_anchor_points = get_genes_hull_points(existing_elements,pos,updated_elements)

    add_signature_metanodes(genes_hull_points,existing_elements,stylesheet_detail,updated_elements,genes_anchor_points)
          


    # print( existing_elements,stylesheet_detail,{"name":"preset",'positions': {
    #         node['data']['id']: node['position']
    #         for node in existing_elements if "source" not in node["data"]
    #     }})

    return updated_elements,stylesheet_detail,{"name":"preset",'positions': {
            node['data']['id']: node['position']
            for node in existing_elements if "source" not in node["data"]
        },
        "animate":False},dict([(i["data"]["id"],i) for i in existing_elements if "source" not in i["data"] ])
def redraw(existing_elements,detail_pos_store,AR,stylesheet_detail):
    pos = dict([(node_id,[detail_pos_store[node_id]["position"]["x"],detail_pos_store[node_id]["position"]["y"]]) for node_id in detail_pos_store])
    bb = [np.inf,np.inf,-np.inf,-np.inf]
    for n in pos:
        if bb[0]>pos[n][0]:
            bb[0]=pos[n][0]
        if bb[1]>pos[n][1]:
            bb[1]=pos[n][1]
        if bb[2]<pos[n][0]:
            bb[2]=pos[n][0]
        if bb[3]<pos[n][1]:
            bb[3]=pos[n][1]
    cur_ar = (bb[2]-bb[0])/(bb[3]-bb[1])
    for n in pos:
        pos[n][0]=pos[n][0]*AR/cur_ar
        detail_pos_store[n]["position"]["x"] = pos[n][0]
    genes_hull_points,genes_anchor_points = get_genes_hull_points(detail_pos_store.values(),pos,existing_elements)
    add_signature_metanodes(genes_hull_points,existing_elements,stylesheet_detail,detail_pos_store,genes_anchor_points,True)
    for elem in existing_elements:
        if elem["data"]["id"] in pos:
            elem["data"]["pos"]=pos[elem["data"]["id"]]
            elem["position"]={"x":pos[elem["data"]["id"]][0],"y":pos[elem["data"]["id"]][1]}
    return list(detail_pos_store.values())+[i for i in existing_elements if "source" in i["data"]],stylesheet_detail,{"name":"preset",'positions': pos,
        "animate":False},detail_pos_store
    