import dash_cytoscape as cyto
import numpy as np
import pandas as pd
import nx_layout
import tlp_layout
from controller import Controller
from scipy.spatial import ConvexHull

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
            "width":"100%",
            "height":"100%","borderWidth":1,"borderColor":"black","borderStyle":"solid"
            },
        autoungrabify=True,
        wheelSensitivity=0.5,
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


def circle_approx_points(center,radius,n_points):
    angles = np.linspace(0,np.pi*2,n_points)
    return np.stack([np.cos(angles)*radius+center[0],np.sin(angles)*radius+center[1]],-1).tolist()

def update_cur_elements(existing_elements,updated,stylesheet_detail,sizes):
    to_del=[]
    for i in range(len(existing_elements)):
        elem = existing_elements[i]
        if not "source" in elem["data"] and elem["data"]["id"] in sizes.index:
            elem["data"]["weight"]=size(sizes.loc[elem["data"]["id"] ]["size"])
            elem["data"]["Signatures"]=sizes.loc[elem["data"]["id"] ]["id"]
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


def add_path_ways(existing_elements,stylesheet_detail):
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

    return existing_elements + list(pathways_nodes.values()) + pathways_edges,stylesheet_detail + list(pathways_stylesheets.values())


def get_genes_hull_points(existing_elements,pos):
    genes_hull_points = {}

    for elem in existing_elements:
        if "source" not in elem["data"]:
            p = pos[elem["data"]['id']]
            if "position" in elem:
                elem["position"]["x"]=p[0]
                elem["position"]["y"]=p[1]
            else:
                elem["position"]={'x':p[0],'y':p[1]}
            elem["data"]["pos"]=elem["position"]
            added_points = circle_approx_points(p,elem["data"]['weight'],12)
            if  "Signatures" in elem["data"]:
                for i in elem["data"]["Signatures"]:
                    if i in genes_hull_points:
                        genes_hull_points[i] = genes_hull_points[i] + added_points
                    else:
                        genes_hull_points[i] = added_points
    return genes_hull_points


def add_signature_metanodes(gene_hull_points,existing_elements,stylesheet_detail):
    for i in gene_hull_points:
        sign_pos = np.array(gene_hull_points[i])
        hull = ConvexHull(sign_pos)
        bb = (np.min(sign_pos[hull.vertices],axis=0),np.max(sign_pos[hull.vertices],axis=0))
        shape_polygon_points  = []

        for v in hull.vertices:
            shape_polygon_points.append(-1 + 2*(sign_pos[v,0]-bb[0][0])/(bb[1][0]-bb[0][0]))
            shape_polygon_points.append(-1 + 2*(sign_pos[v,1]-bb[0][1])/(bb[1][1]-bb[0][1]))
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
        

def color_metanodes(cm,stylesheet_detail):
    for d in cm:
        color_str = f"rgba({','.join([str(i*255) for i in cm[d][:3]])},0.25)"

        stylesheet_detail.append({
            "selector":f"node.{d}",
            "style":{
                "background-color":color_str
            }
        })

def display_detail_graph(selectedDiseases,selected_signatures,selected_genes,existing_elements):
    #TODO add metanodes siwth style shape + label
        
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
    color_by_diseases = True # TODO
    cm = Controller._instance.dm.get_disease_cmap() if color_by_diseases else Controller._instance.dm.get_stage_cmap()

    # print(updated_elements)
    # print(existing_elements)
    # updated_elements.clear()
    print("start")

    intersections,items = Controller._instance.dm.get_genes_intersections(id_filter=selected_signatures,disease_filter=selectedDiseases,gene_filter=selected_genes)

    sizes = items.set_index("gene")
    updated = set()
    
    update_cur_elements(existing_elements,updated,stylesheet_detail,sizes)

        
    nodes = [{'data':{'id':g.gene,"label":g.gene,"weight":size(g.size),"Signatures":g.id}} for g in items.itertuples() if g.gene not in updated]

    edges = [{"data":{"source":k[0],"target":k[1]}} for k,v in intersections.items() ]

    # print(edges)
    existing_elements = existing_elements +nodes
    existing_elements = existing_elements +edges

    existing_elements,stylesheet_detail = add_path_ways(existing_elements,stylesheet_detail)

    pos = tlp_layout.sgd2_layout(existing_elements)

    genes_hull_points = get_genes_hull_points(existing_elements,pos)

    add_signature_metanodes(genes_hull_points,existing_elements,stylesheet_detail)
          

    color_metanodes(cm,stylesheet_detail)

    # print( existing_elements,stylesheet_detail,{"name":"preset",'positions': {
    #         node['data']['id']: node['position']
    #         for node in existing_elements if "source" not in node["data"]
    #     }})
    return existing_elements,stylesheet_detail,{"name":"preset",'positions': {
            node['data']['id']: node['position']
            for node in existing_elements if "source" not in node["data"]
        }}
