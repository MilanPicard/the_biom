import dash_cytoscape as cyto
import pandas as pd
from dash import html
# from dash_extensions.enrich import html

# Initialize Cytoscape with extra layouts
cyto.load_extra_layouts()

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
def overview_graph(dm):
    cyto.load_extra_layouts()

    # intersections,signatures_ids = dm.get_signatures_intersections()
    genes = {}
    edges = dict()
    # for item in signatures.itertuples():
    #     for g in item.Signature.split(";"):
    #         if g in genes:
    #             genes[g].append(item.id)
    #         else:
    #             genes[g] = [item.id]
    # for g in genes:
    #     signs = sorted(genes[g])
    #     for id1 in range(len(signs)):
    #         for id2 in range(id1+1,len(signs)):
    #             edge = Edge(signs[id1],signs[id2])
    #             if(edge in edges):
    #                 edges[edge].append(g)
    #             else:
    #                 edges[edge] = [g]
    return cyto.Cytoscape(
        id="overview_graph",
        layout={"name":"cose","nodeDimensionsIncludeLabels":True,"animate":False},
        clearOnUnhover=True,

        className="overview",
        boxSelectionEnabled=True,
        autoungrabify=True,
        wheelSensitivity=0.2,
        maxZoom=4,
        minZoom=0.25,
        # minZoom=1,
        style={"width":"97%","height":"96%"},
        contextMenu=[
            {
                "availableOn":["canvas"],
                "label":"export as Image",
                "id":"overview_export_image",
                "onClickCustom":"export_image_event_handler"
            },
            {
                "availableOn":["canvas"],
                "label":"export as json",
                "id":f"overview_export_text",
                "onClickCustom":"export_text_event_handler"
            }
        ],
        elements=
            # [{'data':{'id':s.id,"genes":s.Signature,"label":s.id,"Disease":s.Disease}} for s in signatures.itertuples()]+
            # [e.data for e in edges]
            # [{'data':{"id":i["id"],"label":i["id"],"Disease":i["Disease"]}} for i in signatures_ids]+[{"data":{"source":k[0],"target":k[1]}} for k,v in intersections.items()]
            get_elements(dm)
    )
def get_node(elem):
    comp_fields = elem["Comparison"].split("vs")
    label_fields = [elem["Cancer"],comp_fields[1]
                    ]
    return {'data':{"id":elem["id"],"label":"\n".join(label_fields),"Cancer":elem["Cancer"],"Comparison":elem["Comparison"],"Filter":elem["Filter"],"Signature":elem["Signature"],"tooltip_content":[
        html.H6(elem["id"]),
        html.Button("Copy genes to clipboard",id={"type":"signature_clipboard","sign":elem["id"]} ,value=";".join(elem["Signature"]),className="btn btn-info"),
        html.Br(),
        html.A("gProfiler",href=elem["gProfiler"],target="_blank")

    ],"gProfiler":elem["gProfiler"]},"group":"nodes","classes":" ".join([elem["Cancer"]])}
def get_elements(dm,**dm_kwargs):
    intersections,signatures_ids,pathway_edges,_ = dm.get_signatures_intersections(**dm_kwargs)
    genes = set()
    for v in intersections.values():
        genes.update(v)
    symbols = dm.get_symbol(genes)
    cancers = {}
    for i in signatures_ids:
        if i["Cancer"] not in cancers:
            cancers[i["Cancer"]] = [i["id"]]
        else:
            cancers[i["Cancer"]].append(i["id"])

    fake_edges = []
    fake_nodes = []
    for c,l in cancers.items():
        fake_nodes.append({'data':{"id":str(c),"label":"","fake":True},"group":"nodes","style":{"width":0,"height":0}})
        for i in range(len(l)):
            fake_edges.append({"data":{"source":l[i],"target":c,"fake":True,"type":"fake"},"group":"edges","style":{"width":0}})
    
    edges = []
    for k, v in intersections.items():
        edge_symbols = []
        for gene_id in v:
            if gene_id in symbols:
                edge_symbols.extend(symbols[gene_id])
        
        tooltip_string = ", ".join(edge_symbols)
        
        edge_data = {
            "source": k[0],
            "target": k[1],
            "elems": v,
            "symbols": tooltip_string,
            "type": "signature"
        }
        edges.append({"data": edge_data, "group": "edges", "classes": "", "style": {"width": 5 + len(v)}})

    # DEBUG: Print the first edge's data to check the tooltip string
    if edges:
        # print("DEBUG [overview.py]: First edge data created:", edges[0])
        pass

    return [get_node(i) for i in signatures_ids] + edges + fake_nodes + fake_edges + [{
        "data":{
            "source":i.Index.split("__")[0],
            "target":i.Index.split("__")[1],
            "symbols":[", ".join(j) for j in i.PathwayStId],
            "type":"pathway"
        },
        "classes":" ".join([
            "pathway"
        ]),
    } for i in pathway_edges.itertuples()]

def get_default_stylesheet(dm,color_by_diseases=True):
    cm = dm.get_disease_cmap()# if color_by_diseases else dm.get_stage_cmap()
    
    s= [
        {"selector":"node","style":{"label":"data(label)","text-wrap":"wrap","background-opacity":0.25,"width":50,"height":50,}},
        {"selector":"node.highlight","style":{"background-opacity":1}},
        {"selector":"node:selected","style":{"background-opacity":1,"border-color":"black"}},
        {"selector":"node.half_highlight","style":{"background-opacity":0.5}},
        {"selector":"edge","style":{"line-opacity":0.5}},
        {"selector":"edge.highlight","style":{"line-opacity":1}},
        {"selector":"edge.half_highlight","style":{"line-opacity":0.75}},
        {"selector":"edge.pathway","style":{"line-style":"solid","line-dash-pattern":[3,6],"line-opacity":0.25,"line-color":"blue","curve-style":"bezier",  "control-point-distances": 120,  "control-point-weights": 0.1,"control-point-step-size":10,}},
        {"selector":".edge-hover","style":{"line-color":"black", "line-opacity": 0.8, "transition-property":"line-color, line-opacity","transition-duration":"0.2s"}},
    {
        "selector":"node.testFlash",
        "style":{
            "borderColor":"yellow",
            "borderWidth":"2px",
            }
    },
    {
        "selector":"edge.testFlash",
        "style":{
            "lineStyle":"dotted",
            }
    }
            ]
    for d in cm:
        color_str = f"rgba({','.join([str(i*255) for i in cm[d][:3]])},{cm[d][3]})"
        s.append({
            "selector":f"node.{d}",
            "style":{
                "background-color":color_str
            }
        })

    return s
    