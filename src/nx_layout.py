import networkx as nx
import numpy as np
def spring_layout(elements):
    G = nx.Graph()
    min_w = np.min([i["data"]["weight"] for i in elements if "source" not in i["data"]])
    max_w = np.max([i["data"]["weight"] for i in elements if "source" not in i["data"]])
    G.add_nodes_from([i["data"]["id"] for i in elements if "source" not in i["data"]])
    pos = None
    if all(["x" in i["data"] and "y" in i["data"] for i in elements if "source" not in i["data"]]):
        pos = dict([(id["data"]["id"],(i["data"]["x"],i["data"]["y"])) for i in elements if "source" not in i["data"]])
    G.add_edges_from([(i["data"]["source"],i["data"]["target"]) for i in elements if "source" in i["data"]])
    return nx.spring_layout(G,scale=200*max_w/min_w)
