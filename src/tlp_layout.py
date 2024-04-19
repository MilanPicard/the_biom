import numpy as np
from tulip import tlp
import networkx as nx
import s_gd2
def sgd2_layout(elements):
    G = tlp.newGraph()
    min_w = np.min([i["data"]["weight"] for i in elements if "source" not in i["data"]])
    max_w = np.max([i["data"]["weight"] for i in elements if "source" not in i["data"]])
    node_dict = dict()
    pos = {}
    vs = G.getSizeProperty("viewSize")

    for i in elements:
        if "source" not in i["data"]:
            node_dict[i["data"]["id"]]=G.addNode()
            if 'x' in i["data"] and "y" in i["data"]:
                pos[i["data"]["id"]]=(i["data"]["x"],i['data']["y"],0.0)
            vs[node_dict[i["data"]["id"]]] = (4*i["data"]["weight"],4*i["data"]["weight"],1)
        else:
            G.addEdge(node_dict[i["data"]["source"]],node_dict[i["data"]["target"]])
    vl = G.getLayoutProperty("viewLayout")
    ccs = tlp.ConnectedTest.computeConnectedComponents(G)
    for cc in ccs:
        sG = G.inducedSubGraph(cc)
        node_index = dict([(n,i) for i,n in enumerate(cc)])
        I=[]
        J=[]
        for e in sG.edges():
            src = sG.source(e)
            tgt = sG.target(e)
            I.append(node_index[src])
            J.append(node_index[tgt])
        
        X = s_gd2.layout(I,J)
        X*=100*sG.numberOfEdges()/sG.numberOfNodes()
        for n,i in node_index.items():
            vl[n]=(X[i,0],X[i,1],0.0)
    # if len(pos)==len(node_dict):
    #     params["has initial layout"]=True
    #     for i in node_dict:
    #         n = node_dict[i]
    #         vl[n] = pos[i]
            
    params = {'result': vl,
               'termination criterion': 'none',
                'fix x coordinates': False, 
                'fix y coordinates': False, 
                'fix z coordinates': False, 
                'has initial layout': False, 
                'layout components separately': False, 
                'number of iterations': 200, 
                'edge costs': 150.0, 
                'use edge costs property': False, 
                'edge costs property': None}
    # G.applyLayoutAlgorithm("Stress Minimization (OGDF)",params)
    params = {'result': vl, 'stop tolerance': 0.001, 'used layout': len(pos)==len(node_dict), 'zero length': 0.0, 'edge length': 100.0,
               'compute max iterations': True, 'global iterations': 50, 'local iterations': 50}
    # G.applyLayoutAlgorithm("Kamada Kawai (OGDF)",params)
    params = {'result': vl, '3D layout': False, 'edge length': 1000, 'initial layout': vl if len(pos)==len(node_dict) else None, 'unmovable nodes': None, 'max iterations': 0}
    # G.applyLayoutAlgorithm("GEM (Frick)",params=params)
    
    params = {'result': vl, 'initial layout': vl, 'node size': vs, 'rotation': G.getDoubleProperty("viewRotation"), 'complexity': 'auto'}
    G.applyLayoutAlgorithm("Connected Components Packing",params)

    for i in node_dict:
        p = vl[node_dict[i]]
        pos[i]=(p[0],p[1])
    return pos

def force_layout(elements):
    G = tlp.newGraph()
    min_w = np.min([i["data"]["weight"] for i in elements if "source" not in i["data"]])
    max_w = np.max([i["data"]["weight"] for i in elements if "source" not in i["data"]])
    node_dict = dict()
    pos = {}
    vs = G.getSizeProperty("viewSize")

    for i in elements:
        if "source" not in i["data"]:
            node_dict[i["data"]["id"]]=G.addNode()
            if 'x' in i["data"] and "y" in i["data"]:
                pos[i["data"]["id"]]=(i["data"]["x"],i['data']["y"],0.0)
            vs[node_dict[i["data"]["id"]]] = (4*i["data"]["weight"],4*i["data"]["weight"],1)
        else:
            G.addEdge(node_dict[i["data"]["source"]],node_dict[i["data"]["target"]])
    # print(G.numberOfEdges())
    vl = G.getLayoutProperty("viewLayout")

    if len(pos)==len(node_dict):
        params["has initial layout"]=True
        for i in node_dict:
            n = node_dict[i]
            vl[n] = pos[i]
            
    params = {'result': vl,
               'termination criterion': 'none',
                'fix x coordinates': False, 
                'fix y coordinates': False, 
                'fix z coordinates': False, 
                'has initial layout': False, 
                'layout components separately': False, 
                'number of iterations': 200, 
                'edge costs': 150.0, 
                'use edge costs property': False, 
                'edge costs property': None}
    G.applyLayoutAlgorithm("Stress Minimization (OGDF)",params)
    params = {'result': vl, 'stop tolerance': 0.001, 'used layout': len(pos)==len(node_dict), 'zero length': 0.0, 'edge length': 100.0,
               'compute max iterations': True, 'global iterations': 50, 'local iterations': 50}
    # G.applyLayoutAlgorithm("Kamada Kawai (OGDF)",params)
    params = {'result': vl, '3D layout': False, 'edge length': 1000, 'initial layout': vl if len(pos)==len(node_dict) else None, 'unmovable nodes': None, 'max iterations': 0}
    # G.applyLayoutAlgorithm("GEM (Frick)",params=params)
    
    params = {'result': vl, 'initial layout': vl, 'node size': vs, 'rotation': None, 'complexity': 'auto'}
    G.applyLayoutAlgorithm("Connected Components Packing",params)

    for i in node_dict:
        p = vl[node_dict[i]]
        pos[i]=(p[0],p[1])
    return pos
