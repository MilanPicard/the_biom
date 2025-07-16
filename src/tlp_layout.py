import numpy as np
from tulip import tlp
import networkx as nx
import s_gd2
def initialized_sgd2_layout(I,J,V,pos):
    seed = s_gd2.s_gd2._check_random_seed()
    X = s_gd2.s_gd2.random_init(I,J,seed,None)
    for i in pos:
        X[i][0]=pos[i][0]
        X[i][1]=pos[i][1]
    if(len(pos)!=len(X)):
        s_gd2.s_gd2.cpp.layout_weighted(X,I,J,V,50,0.01,seed)
        # X*=100*len(I)/len(X)
        average_edge_length = (((X[I]-X[J])**2).sum(-1)**0.5).mean()
        X*=average_edge_length*50
    return X
def sgd2_layout(elements,detail_pos_store,AR=1):
    G = tlp.newGraph()
    node_dict = dict()
    pos = {}
    vs = G.getSizeProperty("viewSize")
    cur_pos=dict()
    w = G.getDoubleProperty("edgeWeight")
    for i in elements:
        if "source" not in i["data"] and ("is_metanode" not in i["data"] or not i["data"]["is_metanode"]):
            node_dict[i["data"]["id"]]=G.addNode()
            if 'x' in i["data"] and "y" in i["data"]:
                pos[i["data"]["id"]]=(i["data"]["x"],i['data']["y"],0.0)
            if i["data"]["id"] in detail_pos_store:
                if "position" in detail_pos_store[i["data"]["id"]]:
                    cur_pos[ node_dict[i["data"]["id"]]]=(detail_pos_store[i["data"]["id"]]["position"]["x"],detail_pos_store[i["data"]["id"]]["position"]["y"])
            vs[node_dict[i["data"]["id"]]] = (4*i["data"]["weight"],4*i["data"]["weight"],1)
        elif "source" in i["data"]:
            e = G.existEdge(node_dict[i["data"]["source"]],node_dict[i["data"]["target"]],False)
            if(not e.isValid()):
                e = G.addEdge(node_dict[i["data"]["source"]],node_dict[i["data"]["target"]])
                if not "is_pathway_edge" in i["data"] or not i["data"]["is_pathway_edge"]:
                    w[e]=w[e]+1
                else:
                    w[e]=w[e]+3

    vl = G.getLayoutProperty("viewLayout")
    ccs = tlp.ConnectedTest.computeConnectedComponents(G)
    for cc in ccs:
        sG = G.inducedSubGraph(cc)
        node_index = dict([(n,i) for i,n in enumerate(cc)])
        I=[]
        J=[]
        V = []
        for e in sG.edges():
            src = sG.source(e)
            tgt = sG.target(e)
            I.append(node_index[src])
            J.append(node_index[tgt])
            V.append(w[e])
        
        X = initialized_sgd2_layout(I,J,V,dict([(i,cur_pos[n]) for n,i in node_index.items() if n in cur_pos]))
        
        cur_ar = (X[:,0].max()-X[:,0].min())/(X[:,1].max()-X[:,1].min())
        X[:,0]*=AR/cur_ar
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
