import dash_cytoscape as cyto
from dash import Patch
import dash
import numpy as np
from data_manager import DataManager
from scipy.spatial import ConvexHull
import plotly.express as px
from shapely import Polygon,Point
from dash import html,ctx
from dash_svg import Svg, G, Path, Circle,Rect
import euler_layout
from urllib.parse import quote
import urllib.parse
# Module-level variable to store the last green genes
last_green_genes = []

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
def detail_graph(elemId):
    return cyto.Cytoscape(
        id=elemId,
        style={
            "width":"97%",
            "height":"100%"},
        boxSelectionEnabled=True,
        autoungrabify=False,
        wheelSensitivity=0.2,
        layout={"name":"preset"},
        clearOnUnhover=True,
        contextMenu=[
            {
                "availableOn":["canvas"],
                "label":"export as Image",
                "id":f"{elemId}_export_image",
                "onClickCustom":"export_image_event_handler"
            },
            {
                "availableOn":["canvas"],
                "label":"export as json",
                "id":f"{elemId}_export_text",
                "onClickCustom":"export_text_event_handler"
            },{
                "availableOn":["node"],
                "label":"visit gProfiler/Reactome/Ensembl",
                "id":f"{elemId}_node_ctx",
                "onClickCustom":"node_event_handler"
            }
            
        ],
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
    return np.sqrt(2*area)*2

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

def add_path_ways(existing_elements,stylesheet_detail,updated_elements,dm,legend,show_pathways=True,selected_signatures=None):
    if not show_pathways:
        # Do not add any pathway nodes or edges
        return existing_elements, stylesheet_detail, None
    # Decide whether to keep all pathways or only those linking multiple signatures
    pathways_nodes = dict()
    protein_nodes = dict()
    pathways_stylesheets = dict()
    pathways_edges = []
    gene_sign = dict([(elem['data']['id'],set(elem['data']["Signatures"])) for elem in existing_elements if "source" not in elem["data"]])
    # Determine if this is mono view based on the number of signatures being displayed
    if selected_signatures is not None:
        all_pathway = len(selected_signatures) <= 1
    else:
        # Fallback: If only one signature is present among the current genes, show all pathways (mono view)
        unique_signatures = set()
        for signs in gene_sign.values():
            unique_signatures.update(signs)
        all_pathway = len(unique_signatures) <= 1
    pathways = dm.get_pathways(gene_sign.keys())
    pathway_shape = "triangle"
    multi_pathway_color = "red"
    mono_pathway_color = "gray"
    for p in pathways.itertuples():
        children = []
        children.append(html.H6(p.PathwayDisplayName))
        children.append(html.A("PathwayReactomeLink",href=p.PathwayReactomeLink,target="_blank"))
        pathways_nodes[p.Index]={"data":{"id":p.Index,"label":"","weight":len(p.EnsemblID),"tooltip_content":children,"is_pathways":True,"ReactomeLink":p.PathwayReactomeLink,"name":p.PathwayDisplayName},"selectable":True,"classes":"pathway"}
        pathways_stylesheets[p.Index] = {
                    "selector":f"node#{p.Index}",
                    "style":{
                        "shape":pathway_shape,
                        "z-index":-0,
                        "width":size_triangle(len(p.EnsemblID)),
                        "height":size_triangle(len(p.EnsemblID))
                        }
                    }
    for p in pathways.itertuples():
        highlight = np.zeros(len(p.EnsemblID),dtype=bool)
        i=0
        while(i<len(p.EnsemblID)):
            j=0
            while(not highlight[i] and j<len(p.EnsemblID)):
                if i==j:
                    j+=1
                    continue
                sign_i = gene_sign[p.EnsemblID[i]]
                sign_j = gene_sign[p.EnsemblID[j]]
                if(len(sign_i.symmetric_difference(sign_j))!=0):
                    highlight[i]=True
                    highlight[j]=True
                j+=1
            i+=1
        for i in range(len(p.EnsemblID)):
            edge_id = "_".join([p.EnsemblID[i],p.Index])
            pathways_edges.append({"data":{"id":edge_id,"source":p.EnsemblID[i],"target":p.Index,"is_pathway_edge":True,"highlight":1 if highlight[i] else 0},"classes":"pathway"})
    legend["pathway"]={}
    if not all_pathway:
        kept_pathways = set()
        for e in pathways_edges:
            if e["data"]["highlight"]==1:
                kept_pathways.add(e["data"]["target"])
        pathways_edges = [e for e in pathways_edges if e["data"]["target"] in kept_pathways]
        pathways_nodes = {k:v for k,v in pathways_nodes.items() if k in kept_pathways}
        legend["pathway"]["pathways linking multiple signatures"]={"shape":pathway_shape,"color":multi_pathway_color}
        pathways = pathways.loc[pathways.index.isin(kept_pathways)]
    else:
        legend["pathway"]["pathways"]={"shape":pathway_shape,"color":mono_pathway_color}
        # In mono view, keep all pathways and their connections
        # pathways data is already complete, no filtering needed
    return existing_elements + list(pathways_nodes.values())  + list(protein_nodes.values())+ pathways_edges,stylesheet_detail + list(pathways_stylesheets.values()), pathways


def get_genes_hull_points(existing_elements,pos,updated_elements):
    genes_hull_points = {}
    genes_anchor_points = {}

    for i,elem in enumerate(existing_elements):
        if "source" not in elem["data"] and ("is_metanode" not in elem["data"] or not elem["data"]["is_metanode"]):
            p = pos[elem["data"]['id']]
            if "position" in elem:
                elem["position"]["x"]=p[0]
                elem["position"]["y"]=p[1]
                # updated_elements[i]["position"]["x"]=p[0]
                # updated_elements[i]["position"]["y"]=p[1]
            else:
                elem["position"]={'x':p[0],'y':p[1]}
                # updated_elements[i]["position"]={'x':p[0],'y':p[1]}
            elem["data"]["pos"]=elem["position"]
            # updated_elements[i]["data"]["pos"]=elem["position"]
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

def add_signature_metanodes(gene_hull_points,existing_elements,stylesheet_detail,updated_elements,genes_anchor_points,sign_info,style_only=False,convex=True):
    
    index = dict([(elem["selector"],i) for i,elem in enumerate(stylesheet_detail)])

    for i in gene_hull_points:
        sign_pos = np.array(gene_hull_points[i])
        # anchor_pos = np.array(genes_anchor_points[i])
        if convex:
            hull = ConvexHull(sign_pos)
            vertices = hull.vertices
        else:
            vertices = np.arange(len(sign_pos))
        bb = (np.min(sign_pos[vertices],axis=0),np.max(sign_pos[vertices],axis=0))
        shape_polygon_points  = []
        # shape_polygon_anchors  = []
        width = bb[1][0]-bb[0][0]
        height = bb[1][1]-bb[0][1]

        for v in vertices:
            shape_polygon_points.append(-1 + 2*(sign_pos[v,0]-bb[0][0])/(width))
            shape_polygon_points.append(-1 + 2*(sign_pos[v,1]-bb[0][1])/(height))
            # shape_polygon_anchors.append(-1 + 2*(anchor_pos[v//(len(sign_pos)//len(anchor_pos)),0]-bb[0][0])/(bb[1][0]-bb[0][0]))
            # shape_polygon_anchors.append(-1 + 2*(anchor_pos[v//(len(sign_pos)//len(anchor_pos)),1]-bb[0][1])/(bb[1][1]-bb[0][1]))
        label_pos = choose_signature_label_pos(shape_polygon_points)
        stylesheet = {
            "selector":f"node#{i}",
            "style":{
                "width":f"{width}",
                "height":f"{height}",
                "shape":"polygon",
                "backgroundOpacity":0.25,
                "fontSize":10,
                "zIndex":-1,
                "textHalign":label_pos[0],
                "textValign":label_pos[1],
                "borderWidth":1,
                "textWrap":"wrap",
                "shapePolygonPoints":" ".join(list(map(str,shape_polygon_points))),
                # "label":"\n".join(i.split("_"))
                }
                }
        if not style_only:    
            split = i.split("_")
            
            # Generate gProfiler URL dynamically for this signature
            # Get genes for this signature from the data
            signature_genes = []
            try:
                # Extract genes that belong to this signature
                for elem in existing_elements:
                    if "source" not in elem["data"] and not elem["data"].get("is_pathways", False):
                        if i in elem["data"].get("Signatures", []):
                            gene_id = elem["data"]["id"]
                            # Use Ensembl ID directly (gProfiler expects Ensembl IDs)
                            if gene_id.startswith("ENSG"):  # Ensure it's an Ensembl ID
                                signature_genes.append(gene_id)
                
                # Generate dynamic gProfiler URL
                gprofiler_url = generate_gprofiler_url(i, signature_genes)
            except Exception as e:
                #print(f"Error generating gProfiler URL for signature {i}: {e}")
                gprofiler_url = ""
            
            existing_elements.append({
                "data":{"id":i,"label":"","is_metanode":True,"tooltip_content":html.Div([
                    html.H6(i),
                    html.A("gProfiler",href=gprofiler_url,target="_blank"),

                    # html.Button("Copy genes to cliboard",id={"type":"signature_clipboard","sign":i} ,value=";".join(test[i].tolist()))
                    ])},
                "position":{"x":0.5*(bb[0][0]+bb[1][0]),
                            "y":0.5*(bb[0][1]+bb[1][1])},
                "selectable": False,"grabbable":False,
                "classes":" ".join([split[0],split[1].split("vs")[1]])})
            # updated_elements.append({
            #     "data":{"id":i,"label":i,"is_metanode":True,"tooltip_content":i},
            #     "position":{"x":0.5*(bb[0][0]+bb[1][0]),
            #                 "y":0.5*(bb[0][1]+bb[1][1])},
            #     "selectable": False,
            #     "classes":" ".join([j for j in i.split("_") if j!="vs"])})
            stylesheet_detail.append(stylesheet)
        else:
            if not stylesheet["selector"] in index:
                stylesheet_detail.append(stylesheet)
            else:
                stylesheet_detail[index[stylesheet["selector"]]]=stylesheet
            updated_elements[i]["position"]={"x":0.5*(bb[0][0]+bb[1][0]),"y":0.5*(bb[0][1]+bb[1][1])}
        

def color_metanodes(cm,stylesheet_detail,unique_diseases,unique_comparisons):
    signature_legend = {"diseases":{},"comparisons":{}}
    for d in cm:
        color_str = f"rgba({','.join([str(i*255) for i in cm[d][:3]])},0.25)"
        if d in unique_diseases:
        
            stylesheet_detail.append({
                "selector":f"node.{d}",
                "style":{
                    "background-color":color_str,
                    # "background-image":f'data:image/svg+xml;utf8,',
                    # "background-image":f'data:image/svg+xml;utf8,{quote(svg)}',
                    # "background-image":dash.get_asset_url("stripe.svg"),
                    # "background-repeat":"repeat"
                }
            })
            signature_legend["diseases"][d]=color_str
    for i in range(1,5):
        angle = (i-1)*45
        match i:
            case 1:
                roman = "I"
                x1 = "0%"
                x2 = "100%"
                y1 = "0%"
                y2 = "0%"
            case 2:
                roman = "II"
                x1 = "0%"
                x2 = "100%"
                y1 = "0%"
                y2 = "100%"
            case 3:
                roman = "III"
                x1 = "0%"
                x2 = "0%"
                y1 = "0%"
                y2 = "100%"
            case 4:
                x1 = "100%"
                x2 = "0%"
                y1 = "0%"
                y2 = "100%"
                roman = "IV"
        if f"Stage{roman}" in unique_comparisons:
            # Original SVG for main graph hulls
            svg = f'<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE svg><svg version="1.1" baseProfile="full" xmlns="http://www.w3.org/2000/svg" width="40" height="40"> <defs> <linearGradient id="grad"  x1="{x1}" x2="{x2}" y1="{y1}" y2="{y2}"> <stop stop-color="white" stop-opacity="0" offset="0%" /> <stop stop-color="white" stop-opacity="0" offset="43%" /> <stop stop-color="black" stop-opacity="0.25" offset="44%" /> <stop stop-color="black" stop-opacity="0.25" offset="50%" /> <stop stop-color="black" stop-opacity="0.25" offset="56%" /> <stop stop-color="white" stop-opacity="0" offset="57%" /> <stop stop-color="white" stop-opacity="0" offset="100%" /> </linearGradient > </defs> <rect width="100%" height="100%" fill="url(#grad)"/></svg>'
            stylesheet_detail.append({
                "selector":f"node.Stage{roman}",
                "style":{
                    "background-image":f'data:image/svg+xml;utf8,{quote(svg)}',
                    "background-repeat":"repeat"
                }
            })
            # Small, high-contrast SVG for the legend swatch
            legend_svg = "data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='6' height='6'><rect width='6' height='6' fill='white'/><rect x='0' y='0' width='3' height='6' fill='black' fill-opacity='0.25' transform='rotate(45 3 3)'/></svg>"
            signature_legend["comparisons"][f"Normal vs Stage{roman}"]={
                "background-image":f'data:image/svg+xml;utf8,{quote(svg)}',
                "legend-background-image": legend_svg,
                "background-repeat":"repeat"
            }
    return signature_legend

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

def display_detail_graph(selectedDiseases,selected_signatures,selected_genes,existing_elements,detail_pos_store,AR,selected_filter,comparison_filter,show_pathways=True,highlight_gene_ids=None):
    global last_green_genes
    updated_elements = Patch()
    legend_data = {}
    edges = []
    nodes = []
    mono_pathway_color = "gray"
    multi_pathway_color = "red"
    
    # Determine if we are showing common pathways (red) or all (gray)
    show_common_pathways = show_pathways  # In current logic, show_pathways means common pathways

    # In mono-signature mode, pathways (triangles) should be gray; when multiple signatures, shared pathways are red with blue edges
    is_mono_signature_mode = len(selected_signatures) <= 1
    stylesheet_detail = [{
        "selector":"node",
        "style":{
            "width":"data(weight)",
            "height":"data(weight)",
            "label":"data(label)",
            "backgroundColor":"black",
            "zIndex":2
        }
    },
    {
        "selector":"node.testFlash",
        "style":{
            "borderColor":"yellow",
            "borderWidth":"4",
            "background-color":"black"
        }
    },
    {
        "selector":"edge.testFlash",
        "style":{
            "lineStyle":"dotted",
        }
    },
    {
        "selector":"node.pathway",
        "style":{
            "backgroundColor": mono_pathway_color if is_mono_signature_mode else multi_pathway_color,
            "borderColor": mono_pathway_color if is_mono_signature_mode else "black",
            "borderWidth":2,
        }
    },
    {
        "selector":"edge.pathway",
        "style":{
            "lineColor": mono_pathway_color if is_mono_signature_mode else "blue",
            "width": 1 if is_mono_signature_mode else 4,
        }
    },
    {
        "selector":"edge.pathway.tapped",
        "style":{
            "lineColor":"red", 
        }
    },
    {"selector": ".edge-hover", "style": {"width": 8, "line-color": "black", "transition-property": "width, line-color", "transition-duration": "0.2s"}},
    ]
    stylesheet_detail = color_selected_node(stylesheet_detail,selected_genes)
    color_by_diseases = True # TODO
    cm = DataManager._instance.get_disease_cmap()# if color_by_diseases else DataManager.get_instance().get_stage_cmap()

    # print(updated_elements)
    # print(existing_elements)
    # updated_elements.clear()
    intersections,items = DataManager._instance.get_genes_intersections(id_filter=selected_signatures,disease_filter=selectedDiseases,gene_filter= [],selected_filter=selected_filter,comparisons_filter=comparison_filter)
    data = DataManager._instance.get_genes_intersection_data(id_filter=selected_signatures,disease_filter=selectedDiseases,gene_filter= [],selected_filter=selected_filter,comparisons_filter=comparison_filter)
    sign_info = data.filter(["gProfiler","id"]).groupby("id").agg(lambda a: a.iloc[0])["gProfiler"].to_dict()
    unique_signatures = data["id"].unique()
    unique_diseases = np.unique([i.split("_")[0] for i in unique_signatures])
    unique_comparisons = np.unique([i.split("_")[1].split("vs")[1] for i in unique_signatures])

    signature_legend = color_metanodes(cm,stylesheet_detail,unique_diseases,unique_comparisons)
    legend_data.update(signature_legend)
    if(items is None):
        # print("redraw")
        return redraw([],{},AR,stylesheet_detail)
    # sizes = items.set_index("gene")
    updated = set()
    # update_cur_elements(existing_elements,updated,stylesheet_detail,sizes,updated_elements)
    symbols = DataManager._instance.get_symbol([g.gene for g in items.itertuples() if g.gene not in updated])
    all_nodes = []
    for g in items.itertuples():
        node_classes = []
        # Mark as green node if in multiple signatures
        if isinstance(g.id, (list, set)) and len(g.id) > 1:
            node_classes.append('multi-signature-gene')
        # Add highlight class if needed
        if highlight_gene_ids is not None and g.gene in highlight_gene_ids:
            node_classes.append('legend-highlight')
        all_nodes.append({
            'data': {
                'id': g.gene,
                'label': "",
                'weight': size(g.size),
                'Signatures': g.id,
                'tooltip_content': [
                    html.H6(symbols[g.gene]),
                    html.A("Ensembl", href=f"https://ensembl.org/Homo_sapiens/Search/Results?q={g.gene};site=ensembl_all;facet_species=Human;facet_feature_type=Gene", target="_blank")
                ]
            },
            'grabbable': False,
            'classes': ' '.join(node_classes) if node_classes else None
        })
    # Print and store green nodes (genes in multiple signatures) immediately after graph generation
    green_genes = [n['data']['id'] for n in all_nodes if n.get('classes') and 'multi-signature-gene' in n['classes']]
    last_green_genes = green_genes
    #print("=== GREEN GENES (genes in multiple signatures, upon graph generation) ===")
    #print(green_genes)
    #print("=== END GREEN GENES ===")
    # nodes = [{'data':{'id':g.gene,"label":"","weight":size(g.size),"Signatures":g.id,"tooltip_content":symbols[g.gene]}} for g in items.itertuples() if g.gene not in updated]

    # edges = [{"data":{"source":k[0],"target":k[1]}} for k,v in intersections.items() ]

    # print(edges)
    # existing_elements = existing_elements +nodes
    # for n in nodes:
    #     updated_elements.append(n)
    # existing_elements= existing_elements +edges
    # for e in edges:
    #    updated_elements.append(e)
    all_nodes_and_edges,stylesheet_detail,pathways = add_path_ways(all_nodes,stylesheet_detail,updated_elements,dm=DataManager._instance,legend=legend_data,show_pathways=show_pathways,selected_signatures=selected_signatures)
    pos,hull = euler_layout.euler_layout(data,all_nodes_and_edges,data.filter(["EnsemblID","id"]).groupby("EnsemblID").agg(lambda a: len(a))["id"].max())
    exportable_data = get_exportable_data(items, pathways)

    # if(len(l)!=prev_len):
    #     prev_len=len(l)
    #     print(len(l))


    # existing_elements,stylesheet_detail = add_path_ways(existing_elements,stylesheet_detail,updated_elements,DataManager.get_instance())
    tmp_nodes = dict()
    tmp_edges = []
    for n in all_nodes_and_edges:
        if "source" not in n["data"] and not ("is_pathways" in n["data"] and n["data"]["is_pathways"]):
            for i in n["data"]["Signatures"]:
                if i not in tmp_nodes:
                    tmp_nodes[i]={'data':{'id':i,"weight":1}}
                else:
                    tmp_nodes[i]['data']["weight"]+=1
                tmp_edges.append({'data':{'source':n['data']['id'],"target":i,"weight":1}})
    
    # pos = tlp_layout.sgd2_layout(all_nodes_and_edges+list(tmp_nodes.values())+tmp_edges,detail_pos_store,AR)
    # print("len(pos)",len(pos))

    # genes_hull_points,genes_anchor_points = get_genes_hull_points(all_nodes_and_edges,pos,updated_elements)
    # print(genes_hull_points.keys())
    # add_signature_metanodes(genes_hull_points,all_nodes_and_edges,stylesheet_detail,updated_elements,genes_anchor_points)
    add_signature_metanodes(hull,all_nodes_and_edges,stylesheet_detail,updated_elements,None,convex=False,sign_info=sign_info)
          
    existing_nodes = {n["data"]["id"]:n for n in all_nodes_and_edges if "source" not in n["data"]}
    to_append = {n["data"]["id"]:n for n in all_nodes_and_edges if "source" not in n["data"]}
    existing_edges = {(n["data"]["source"],n["data"]["target"]):n for n in all_nodes_and_edges if "source"  in n["data"]}

    to_del = []
    def update_elem(elem,existing):
        if isinstance(existing,dict):
            for k in existing:
                update_elem(elem[k],existing[k])
        else:
            elem = existing
    for i,n in enumerate(existing_elements):
        if "source" in n["data"]:
            if (n["data"]["source"],n["data"]["target"]) not in existing_edges:
                to_del.append(i)
            else:
                if "highlight" in n["data"] and n["data"]["highlight"]!=existing_edges[(n["data"]["source"],n["data"]["target"])]["data"]["highlight"]:
                    updated_elements[i]["data"]["highlight"]=existing_edges[(n["data"]["source"],n["data"]["target"])]["data"]["highlight"]
                    updated_elements[i]["classes"]=existing_edges[(n["data"]["source"],n["data"]["target"])]["classes"]
                del existing_edges[(n["data"]["source"],n["data"]["target"])]
        else:
            if(n["data"]["id"] in existing_nodes):
                existing = existing_nodes[n["data"]["id"]]
                for field in existing:
                    # updated_elements[i][field]=existing[field]
                    update_elem(updated_elements[i][field],existing[field])
                updated_elements[i]["selected"]=False
                del to_append[n["data"]["id"]]
            else:
                to_del.append(i)
    for n in reversed(to_del):
        del updated_elements[n]
    for n in to_append.values():
        updated_elements.append(n)
    for n in existing_edges.values():
        updated_elements.append(n)
    multi_sign_gene_color = "limegreen"
    mono_sign_gene_color = "black"
    for n in all_nodes:
        if "weight" not in n["data"]:
            continue  # Skip nodes without weight (e.g., if pathways are hidden)
        if n["data"]["id"] not in selected_genes:
            stylesheet_detail.append({
                "selector":"node#"+n["data"]["id"],
                "style":{
                    "height":n["data"]["weight"],
                    "width":n["data"]["weight"],
                    "background-color":multi_sign_gene_color if n["data"]["weight"]>20 else mono_sign_gene_color,
                    "border-color":"black",
                    "border-width":2
                }
            })
        else:
            stylesheet_detail.append({
                "selector":"node#"+n["data"]["id"],
                "style":{
                    "height":n["data"]["weight"],
                    "width":n["data"]["weight"],
                    "border-color":"black",
                    "border-width":2
                }
            })
    legend_data["genes"] = {"genes": mono_sign_gene_color, "genes in multiple signatures": multi_sign_gene_color}
    if not show_pathways:
        legend_data["genes"]["genes in multiple signatures"]=multi_sign_gene_color
    # print( existing_elements,stylesheet_detail,{"name":"preset",'positions': {
    #         node['data']['id']: node['position']
    #         for node in existing_elements if "source" not in node["data"]
    #     }})
    return updated_elements,stylesheet_detail,{"name":"preset",'positions': {
            node['data']['id']: node['position']
            for node in all_nodes_and_edges if "source" not in node["data"]
        },
        "animate":False},existing_nodes,legend_data,exportable_data#dict([(i["data"]["id"],i) for i in existing_elements if "source" not in i["data"] ])

def get_exportable_data(items, pathways):
    if pathways is None:
        # No pathways to export
        return {}
    exportable_data = {"gene_groups":{},"pathways":{}}
    zones = {}
    for g in items.itertuples():
        signs = g.id
        signs_id = ";".join(signs)
        if signs_id not in zones:
            zones[signs_id]={"name":f"{signs[0]} only" if len(signs)==1 else f"intersection of {','.join(signs)}","genes":[g.gene]}
        else:
            zones[signs_id]["genes"].append(g.gene)
    for z in zones.values():
        exportable_data["gene_groups"][z["name"]]=z["genes"]
    for p in pathways.itertuples():
        if p.Index not in exportable_data["pathways"]:
            exportable_data["pathways"][p.Index]={"name":p.PathwayDisplayName,"Reactome":p.PathwayReactomeLink,"connected_genes":{}}
        for g in p.EnsemblID:
            exportable_data["pathways"][p.Index]["connected_genes"][g]=items.loc[items["gene"]==g]["id"].iloc[0]
    return exportable_data
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
    # genes_hull_points,genes_anchor_points = get_genes_hull_points(detail_pos_store.values(),pos,existing_elements)
    # for k in list(genes_hull_points.keys()):
    #     if k not in detail_pos_store:
    #         del genes_hull_points[k]
    #         del genes_anchor_points[k]
    # add_signature_metanodes(genes_hull_points,existing_elements,stylesheet_detail,detail_pos_store,genes_anchor_points,style_only=True)
    for elem in existing_elements:
        if elem["data"]["id"] in pos:
            elem["data"]["pos"]=pos[elem["data"]["id"]]
            elem["position"]={"x":pos[elem["data"]["id"]][0],"y":pos[elem["data"]["id"]][1]}
    for i in stylesheet_detail:
        if "shapePolygonPoints" in i["style"] : 
            i["style"]["width"]=float(i["style"]["width"])*AR/cur_ar
    return list(detail_pos_store.values())+[i for i in existing_elements if "source" in i["data"]],stylesheet_detail,{"name":"preset",'positions': pos,
        "animate":False},detail_pos_store,dash.no_update,dash.no_update
    
def generate_gprofiler_url(signature_id, genes):
    """
    Generate a gProfiler URL dynamically for a given signature and its genes.
    
    Args:
        signature_id (str): The signature ID (e.g., "LIHC_NormalvsStageI_Merge")
        genes (list): List of Ensembl IDs (e.g., ["ENSG00000128710", "ENSG00000104938"])
    
    Returns:
        str: Complete gProfiler URL
    """
    if not genes:
        return ""
    
    # Use Ensembl IDs directly (space-separated as expected by gProfiler)
    genes_param = " ".join(genes)
    
    # Create gProfiler URL with comprehensive parameters
    gprofiler_url = (
        f"https://biit.cs.ut.ee/gprofiler/gost?"
        f"query={genes_param}&"
        f"organism=hsapiens&"
        f"sources=GO:BP,GO:MF,GO:CC,KEGG,REAC,WP&"
        f"user_threshold=0.05&"
        f"all_results=false&"
        f"ordered=false&"
        f"significant=true&"
        f"no_iea=false&"
        f"domain_scope=annotated&"
        f"measure_underrepresentation=false&"
        f"evcodes=false&"
        f"as_ranges=false&"
        f"background=0&"
        f"domain_size_type=known&"
        f"term_size_filter_min=3&"
        f"term_size_filter_max=500&"
        f"numeric_namespace=ENTREZGENE_ACC&"
        f"pictograms=false&"
        f"min_set_size=3&"
        f"max_set_size=500"
    )
    
    return gprofiler_url
    
# Utility function to retrieve green nodes (genes in multiple signatures) from elements

def get_multi_signature_genes(elements):
    """
    Given a list of elements (nodes), return the list of gene IDs that are colored green (i.e., genes in multiple signatures),
    identified by the 'multi-signature-gene' class.
    """
    green_genes = []
    for elem in elements:
        data = elem.get('data', {})
        classes = elem.get('classes', '')
        if (
            'id' in data and data['id'].startswith('ENSG')
            and 'multi-signature-gene' in (classes or '')
        ):
            green_genes.append(data['id'])
    return green_genes
    
def get_last_green_genes():
    """
    Return the most recently generated list of green node gene IDs (genes in multiple signatures).
    """
    global last_green_genes
    return last_green_genes
    
def update_node_highlight_styles(elements, highlight_gene_ids):
    """
    Given a list of elements and a list of gene IDs, add the 'legend-highlight' class to those nodes.
    Remove the class from all others.
    """
    for elem in elements:
        data = elem.get('data', {})
        if 'id' in data and data['id'] in highlight_gene_ids:
            classes = elem.get('classes', '') or ''
            if 'legend-highlight' not in classes:
                elem['classes'] = (classes + ' legend-highlight').strip()
        else:
            classes = elem.get('classes', '') or ''
            if 'legend-highlight' in classes:
                elem['classes'] = ' '.join([c for c in classes.split() if c != 'legend-highlight'])
    return elements

# Add a style for 'legend-highlight' to be used in the Cytoscape stylesheet
legend_highlight_stylesheet = {
    "selector": ".legend-highlight",
    "style": {
        "border-width": 3,
        "border-color": "#FF8800",  # Subtle orange, can be changed
        "border-opacity": 1,
    }
}
    