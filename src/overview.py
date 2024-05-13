import dash_cytoscape as cyto
import pandas as pd

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

    intersections,signatures_ids = dm.get_signatures_intersections()
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
        layout={"name":"cose","nodeDimensionsIncludeLabels":True},
        className="overview",
        boxSelectionEnabled=True,
        autoungrabify=True,
                wheelSensitivity=0.5,
                maxZoom=4,
                minZoom=0.25,
        # minZoom=1,
        style={"width":"97%","height":"96%"},
        elements=
            # [{'data':{'id':s.id,"genes":s.Signature,"label":s.id,"Disease":s.Disease}} for s in signatures.itertuples()]+
            # [e.data for e in edges]
            # [{'data':{"id":i["id"],"label":i["id"],"Disease":i["Disease"]}} for i in signatures_ids]+[{"data":{"source":k[0],"target":k[1]}} for k,v in intersections.items()]
            get_elements(dm)
    )

def get_elements(dm,**dm_kwargs):
    intersections,signatures_ids = dm.get_signatures_intersections(**dm_kwargs)
    return [{'data':{"id":i["id"],"label":"\n".join(i["id"].split("_")),"Disease":i["Disease"],"Comparison":"_".join(i["Comparison"]),"Signature":i["Signature"]},"group":"nodes","classes":" ".join([i["Disease"],"highlight"])} for i in signatures_ids]+[{"data":{"source":k[0],"target":k[1],"elems":v},"group":"edges","classes":""} for k,v in intersections.items()]

def get_default_stylesheet(dm,color_by_diseases=True):
    cm = dm.get_disease_cmap() if color_by_diseases else dm.get_stage_cmap()
    
    s= [
        {"selector":"node","style":{"label":"data(label)","text-wrap":"wrap","background-opacity":0.25}},
        {"selector":"node.highlight","style":{"background-opacity":1}},
        {"selector":"node.half_highlight","style":{"background-opacity":0.5}},
        {"selector":"edge","style":{"line-opacity":0.25}},
        {"selector":"edge.highlight","style":{"line-opacity":1}},
        {"selector":"edge.half_highlight","style":{"line-opacity":0.5}},

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
    