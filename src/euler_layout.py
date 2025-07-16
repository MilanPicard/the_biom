import time
from tulip import tlp
import pandas as pd
import numpy as np
import networkx as nx
import shapely
import cProfile,pstats,io
import impred
"""
Paolo Simonetto, David Auber, Daniel Archambault. Fully Automatic Visualisation of Overlapping
Sets. Computer Graphics Forum, 2009, 28 (3), pp.967-974. ￿hal-00407269
"""


def get_planar_intersection_graph(grouped_data):
    zones = grouped_data["id"].unique().tolist()
    signatures = set()
    for z in zones:
        signatures.update(z.split(";"))
    signatures = list(signatures)
    value_counts = grouped_data["id"].value_counts().reset_index(drop=False)
    value_counts["size"]=value_counts["id"].str.count(";")+1
    value_counts["count"]=value_counts["size"]*value_counts["count"]
    value_counts = value_counts.filter(["id","count"]).set_index("id")["count"].to_dict()
    # print(grouped_data["id"].value_counts())
    # grouped_data["id"].value_counts().reset_index(drop=False).to_parquet("/tmp/value_counts.pq")
    # ccs = np.array(list(range(len(zones))))
    ccs2 = np.arange(len(zones)).reshape(1,-1).repeat(len(signatures),0)
    # np.save("/tmp/zones.npy",np.array(zones))
    # np.save("/tmp/signatures.npy",np.array(signatures))

    for s in range(len(signatures)):
        for j in range(len(zones)):
            if(signatures[s] not in zones[j].split(";")):
               ccs2[s,j]=-1
    # np.save("/tmp/ccs2.npy",ccs2)
    # cc_classes = [set(zones[i].split(";")) for i in ccs]
    # c = np.zeros((len(ccs),len(ccs)))
    # c[0]=ccs[:,0]!=ccs[:,:]
    # for i in range(len(ccs)):
    #     for j in range(i+1,len(ccs)):
    #         c[i,j]=len(cc_classes[i].intersection(cc_classes[j]))#-p1*u-p2*v
    c2 = np.logical_and(np.logical_and(ccs2[:,np.newaxis]>=0 ,ccs2[:,:,np.newaxis]>=0),ccs2[:,np.newaxis]!=ccs2[:,:,np.newaxis]).astype(int).sum(0)
    I = np.zeros(c2.shape+(2,len(signatures)),dtype=bool)
    for i,z in enumerate(zones):
        for s in z.split(";"):
            j = signatures.index(s)
            I[i,:,0,j]=True
            I[:,i,1,j]=True
    
    size = I.astype(int).sum(-1)
    min_size = size.min(-1)
    max_size = size.max(-1)
    intersection = I.all(-2).astype(int).sum(-1)
    U = min_size-intersection
    V = max_size-intersection-1
    discarded_edges =np.zeros((len(c2),len(c2)),dtype=bool)
    g = tlp.newGraph()
    vla = g.getStringProperty("viewLabel")
    nodes = g.addNodes(len(zones))
    for i,n in enumerate(nodes):
        vla[n]=zones[i]
    # signatures_subgraphs = {s:}
    i,j = np.unravel_index(np.argmax(c2-5e-3*U-1e-3*V),c2.shape)
    while(c2[i,j]>0):
        e = g.addEdge(nodes[i],nodes[j])
        is_planar = tlp.PlanarityTest.isPlanar(g)
        if not is_planar:
            t1 = time.perf_counter()
            g.delEdge(e)
            discarded_edges[i,j]=True
            discarded_edges[j,i]=True
            # c[i,j]=0
            c2[i,j]=0
            c2[j,i]=0
        else:
            # kept_index,to_merge_index = (i,j) if ccs[i]<ccs[j] else (j,i)
            # cc_classes[kept_index]=cc_classes[kept_index].union(cc_classes[to_merge_index])
            # # c_copy=np.copy(c)
            # c[to_merge_index,kept_index]=0
            # # c_copy=np.copy(c)

            # c[kept_index,to_merge_index]=0
            # # c_copy=np.copy(c)

            # # for k1 in range(kept_index):
            # #     if not discarded_edges[k1,kept_index]:
            # #         c[k1,kept_index]=len(cc_classes[k1].intersection(cc_classes[kept_index]))
            # # for k1 in range(kept_index+1,len(ccs)):
            # #     if not discarded_edges[kept_index,k1]:
            # #         c[kept_index,k1]=len(cc_classes[k1].intersection(cc_classes[kept_index]))

            # cc_classes[to_merge_index]=set()
            # c[np.unravel_index((np.argwhere(ccs==ccs[kept_index]).reshape(-1,1)*c.shape[0]+np.argwhere(ccs==ccs[to_merge_index]).reshape(1,-1)).flatten(),c.shape)]=0

            # c[np.unravel_index((np.argwhere(ccs==ccs[to_merge_index]).reshape(-1,1)*c.shape[0]+np.argwhere(ccs==ccs[kept_index]).reshape(1,-1)).flatten(),c.shape)]=0
            # ccs[ccs==ccs[to_merge_index]]=ccs[kept_index]
            for k in range(len(ccs2)):
                if ccs2[k,i]>=0 and ccs2[k,j]>=0:
                    ccs2[k,ccs2[k]==ccs2[k,j]]=ccs2[k,i]
            c2 = np.logical_and(np.logical_and(ccs2[:,np.newaxis]>=0 ,ccs2[:,:,np.newaxis]>=0),ccs2[:,np.newaxis]!=ccs2[:,:,np.newaxis]).astype(int).sum(0)*(1-discarded_edges.astype(int))
            
        # i,j = np.unravel_index(np.argmax(c2),c2.shape)
        i,j = np.unravel_index(np.argmax(c2-5e-3*U-1e-3*V),c2.shape)
    # for s in signatures:
    #     sg = g.inducedSubGraph([n for i,n in enumerate(nodes) if s in zones[i].split(";") ])
    #     assert tlp.ConnectedTest.isConnected(sg)
    #     g.delAllSubGraphs(sg)
    # print(times2,times2/times2.sum())
    # print(times[1:]-times[:-1],(times[1:]-times[:-1])/(times[-1]-times[0]))

    return g,zones,value_counts
def apply_bertault(g,surrounding_edges,n_iter=20,edge_length=0,fixed_nodes=set(),gamma=2.0,optimalDist=5.0,bends=False,flexible_edges=False,debug=False):
    # edge_length=0
    # n=0
    # if g.existProperty("fixed"):
    #     fixed = g.getBooleanProperty("fixed")
    #     vl = g.getLayoutProperty("viewLayout")
    #     for e in g.getEdges():
    #         src = g.source(e)
    #         tgt = g.target(e)
    #         if fixed[src] and fixed[tgt]:
    #             edge_length += vl[src].dist(vl[tgt])
    #             n+=1

    # else:
        # edge_length=0
        # n=1
    # params = tlp.getDefaultPluginParameters('Bertault (OGDF)')
    # params["result"]=g.getLayoutProperty("viewLayout")
    # params["number of iterations"]=n_iter
    # params["impred"]=True
    # params["edge length"]=edge_length
    # g.applyLayoutAlgorithm('Bertault (OGDF)',params)
    # import impred
    # i = impred.Impred(g,3,5,n_iter,None)
    # pos = i.layout()
    # import impred_ocotillo
    # impred_ocotillo.impred_ocotillo(g,surrounding_edges,fixed_nodes,gamma=delta,nIter=n_iter,optimalDist=optimalDist,use_bends=bends,flexible_edges=flexible_edges)
    times = [time.perf_counter()]
    
    i = impred.Impred(g,optimalDist,gamma,n_iter,None,surrounding_edges if surrounding_edges is not None else {},fixed_nodes)
    times.append(time.perf_counter())
    i.layout(debug)
    # times.append(time.perf_counter())
    # times = np.array(times)
    # print("b", times[1:]-times[:-1],(times[1:]-times[:-1])/(times[-1]-times[0]))

def tlp2nx(g):
    G = nx.DiGraph()
    G.add_nodes_from(g.nodes())
    G.add_edges_from([(g.source(e),g.target(e)) for e in g.edges()])
    return G

def apply_planar_layout(g,scale=1):
    # params = tlp.getDefaultPluginParameters('Tree Radial')
    vl = g.getLayoutProperty("viewLayout")
    vs = g.getSizeProperty("viewSize")
    vs.setAllNodeValue(tlp.Vec3f(0.1,0.1,0.1))
    cc= tlp.ConnectedTest.computeConnectedComponents(g)
    if True or len(cc)>1:
        avg=[]
        sizes=[]
        for i in range(len(cc)):
            comp = g.inducedSubGraph(cc[i])
            G = tlp2nx(comp)
            pos = nx.planar_layout(G,scale=scale)
            for n,p in pos.items():
                vl[n]=tlp.Vec3f(p[0],p[1],0.0) 
                # print(vs[n])
            if comp.numberOfNodes()>1:
                # if comp.numberOfEdges()>1:
                #     s = 1/vl.averageEdgeLength(comp)
                #     vl.scale(tlp.Vec3f(s,s,1.0),comp)
                delta = vl.averageEdgeLength() if g.numberOfEdges()>0 else 5
                vl.center(comp)
                # c,r = tlp.computeBoundingRadius(comp)
                # delta *=3
                # delta=2
                # vl.scale(tlp.Vec3f(np.sqrt(comp.numberOfNodes())*delta/r.dist(c),np.sqrt(comp.numberOfNodes())*delta/r.dist(c),1.0))

                gamma = delta/2
                edges = comp.edges()
                apply_bertault(comp,{n : edges for n in comp.nodes()} ,n_iter=30,optimalDist=delta,gamma=gamma,bends=False)
            if(comp.numberOfEdges()>1):
                sizes.append(len(cc[i]))
                avg.append(vl.averageEdgeLength(comp))
        if len(sizes)>0:
            sizes = np.array(sizes)
            avg = (np.array(avg)*sizes).sum()/np.sum(sizes)
        g.applyLayoutAlgorithm('Connected Components Packing',{'result': vl, 'initial layout': vl, 'node size': None, 'rotation': None, 'complexity': 'auto'})
        if g.numberOfEdges()>1:
            s = avg/vl.averageEdgeLength(g)
            vl.scale(tlp.Vec3f(s,s,1.0))
    else:
        G = tlp2nx(g)
        pos = nx.planar_layout(G,scale=scale)
        for n,p in pos.items():
            vl[n]=tlp.Vec3f(p[0],p[1],0.0) 
    # params["result"]=vl
    # g.applyLayoutAlgorithm('Random layout')
    # g.applyLayoutAlgorithm('Tree Radial',params)
    vl.center()
def build_grid_graph(g):
    return None
def point_segment_distance(points,srcs,dsts):
    v = dsts-srcs
    d = np.linalg.norm(v,axis=-1)
    unit = v/np.expand_dims(d,-1)
    dot = np.sum( (points-srcs)* unit,-1,keepdims=True)
    r= unit * dot+srcs
    on_edge = np.logical_and(np.logical_or(np.abs(srcs[:,0]-dsts[:,0])<1e-4,np.logical_xor(srcs[:,0]<=r[:,0],dsts[:,0]<=r[:,0])),
                                    np.logical_or(np.abs(srcs[:,1]-dsts[:,1])<1e-4,np.logical_xor(srcs[:,1]<=r[:,1],dsts[:,1]<=r[:,1])))
    
    return np.where(on_edge,np.linalg.norm(points-r,axis=1),np.minimum(np.linalg.norm(points-srcs,axis=1),np.linalg.norm(points-dsts,axis=1)))
def compute_circle_set_radius(g):
    vl = g.getLayoutProperty("viewLayout")
    pos = []
    for n in g.getNodes():
        p = vl[n]
        pos.append((p[0],p[1]))
    pos = np.array(pos)
    if len(pos)>1:
        sd = ((np.expand_dims(pos,0)-np.expand_dims(pos,1))**2).sum(-1)
        np.fill_diagonal(sd,np.inf)
        R = (np.min(sd)**0.5)/3
    else:
        R=10
    if g.numberOfEdges()>0:
        points = []
        srcs = []
        dsts = []
        for e in g.getEdges():
            src = g.source(e)
            tgt = g.target(e)
            mask = np.ones(len(pos),dtype=bool)
            mask[src.id]=False
            mask[tgt.id]=False
            p = vl[src]
            srcs.append(np.array([[p[0],p[1]]]).repeat(len(pos)-2,0))
            p = vl[tgt]
            dsts.append(np.array([[p[0],p[1]]]).repeat(len(pos)-2,0))
            points.append(pos[mask])
        R = min(R,np.min(point_segment_distance(np.concatenate(points,0),np.concatenate(srcs),np.concatenate(dsts))/3))
    return R
def compute_angle(pos_n,pos_neighbor):
    v = pos_neighbor-pos_n
    v = v/np.linalg.norm(v)
    # c = np.dot(v,np.array((1,0)))
    # a = np.arccos(c) +(np.pi if v[1]<0 else 0)
    
    a = np.arctan2(v[1],v[0])
    return a
def construct_hull(zone,zones,g,grid_graph,adjacent_edges_crossing,grid_nodes):
    nodes = g.nodes()
    zone_prop = grid_graph.getStringVectorProperty("zones")
    subgraph_nodes = []
    grid_subgraph_nodes = []
    point = None
    for i,z in enumerate(zones):
        if zone in z.split(";"):
            subgraph_nodes.append(nodes[i])
            grid_subgraph_nodes+=grid_nodes[i]
            for n in grid_nodes[i]:
                zone_prop[n] = sorted(list(set(zone_prop[n]+z.split(";"))))
            if(zone == z):
                point=(grid_nodes[i][len(grid_nodes[i])//2])
    if point is None:
        for i,z in enumerate(zones):
            if zone in z.split(";"):
                point=(grid_nodes[i][len(grid_nodes[i])//2])
                break
    g_subgraph = g.inducedSubGraph(subgraph_nodes,name=zone)
    grid_subgraph = grid_graph.inducedSubGraph(grid_subgraph_nodes,name=f"zone_{zone}")
    vl = grid_subgraph.getStringProperty("viewLabel")
    # vl[point]=zone
 

    for e in g_subgraph.edges():
        if e in adjacent_edges_crossing:
            src = g_subgraph.source(e)
            tgt = g_subgraph.target(e)
            e2 = grid_subgraph.existEdge(*adjacent_edges_crossing[e][src])
            assert e2.isValid()
            grid_subgraph.delEdge(e2)

            e2 = grid_subgraph.existEdge(*adjacent_edges_crossing[e][tgt])
            assert e2.isValid()
            grid_subgraph.delEdge(e2)
    return set(grid_subgraph.edges())

def construct_regions(g,radius,grid_graph,zones,value_counts):
    alpha= 2*np.pi/3
    vl = g.getLayoutProperty("viewLayout")
    vlab = g.getStringProperty("viewLabel")
    vlab_grid = grid_graph.getStringProperty("viewLabel")
    vl_grid = grid_graph.getLayoutProperty("viewLayout")
    vl.computeEmbedding()
    positions = []
    nodes = g.nodes()
    node_sectors = []
    adjacent_edges_crossing={}
    grid_nodes = []
    min_value_count = np.min(list(value_counts.values()))
    for i,n in enumerate(nodes):
        sectors = []
        pos_n = vl[n]
        pos_n = np.array((pos_n.x(),pos_n.y()))
        edges = g.allEdges(n)
        for edge in edges:
            neighbor = g.opposite(edge,n)
            pos_neighbor = vl[neighbor]
            pos_neighbor = np.array((pos_neighbor.x(),pos_neighbor.y()))
            angle = compute_angle(pos_n,pos_neighbor)
            while len(sectors)>0 and angle < sectors[-1]:
                angle += 2*np.pi
                # assert angle >sectors[-1],f"{angle} {sectors}"
            sectors.append(angle)
        if len(sectors)==0:
            sectors.append(0)
        assert sectors[-1]>=sectors[0] , f"{sectors}"
        assert sectors[-1]-sectors[0]<=2*np.pi , f"{sectors}"
        sectors.append(np.pi*2+sectors[0])
        angles = []
        # alpha2 = alpha/(np.floor(1+np.log2((1+value_counts[zones[i]])/(1+1))))
        alpha2 = alpha/(np.floor(1+np.log((1+value_counts[zones[i]])/(1+1))/np.log(4)))
        for i in range(len(sectors)-1):
            start = sectors[i]
            end = sectors[i+1]
            if end<start:
                start-=np.pi*2
            central_angle = end-start
            p = np.ceil(central_angle/alpha2).astype(int)
            angles.append(np.linspace(start,end,2*p+1)[1::2])
        new_nodes = []
        # new_nodes3 = []
        for i in range(len(angles)):
            new_nodes.append(grid_graph.addNodes(len(angles[i])))
            # new_nodes3.append(g.addNodes(len(angles[i])))
            pos = np.expand_dims(pos_n,0)+np.stack((np.cos(angles[i]),np.sin(angles[i])),-1)*radius
            for j in range(pos.shape[0]):
                vl_grid[new_nodes[i][j]]=tlp.Vec3f(pos[j,0],pos[j,1],0.0)
                # vl[new_nodes3[i][j]]=tlp.Vec3f(pos[j,0],pos[j,1],0.0)
                # vlab[new_nodes3[i][j]]=str(i)
        for i,e in enumerate(edges):
            if e not in adjacent_edges_crossing:
                adjacent_edges_crossing[e]={}
            adjacent_edges_crossing[e][n]=(new_nodes[i-1][-1],new_nodes[i][0])
        new_nodes2 = []
        # new_nodes4 = []
        for n2 in new_nodes:
            for n1 in n2:
                new_nodes2.append(n1)
        grid_nodes.append(new_nodes2)
        # for n2 in new_nodes3:
            # for n1 in n2:
                # new_nodes4.append(n1)
        for i in range(len(new_nodes2)):
            grid_graph.addEdge(new_nodes2[i],new_nodes2[(i+1)%len(new_nodes2)])
            vlab_grid[new_nodes2[i]]=vlab[n]
        # for i in range(len(new_nodes4)):
            # g.addEdge(new_nodes4[i],new_nodes4[(i+1)%len(new_nodes4)])
        
    for e in g.getEdges():
        src = g.source(e)    
        tgt = g.target(e)
    #     # if e in adjacent_edges_crossing:
        s1,s2 = adjacent_edges_crossing[e][src]
        t1,t2 = adjacent_edges_crossing[e][tgt]
    #     assert grid_graph.hasEdge(s1,s2)
    #     assert grid_graph.hasEdge(t1,t2)
    #     # grid_graph.delEdge(grid_graph.existEdge(s1,s2))
    #     # grid_graph.delEdge(grid_graph.existEdge(t1,t2))
        grid_graph.addEdge(s1,t2)
        grid_graph.addEdge(t1,s2)
    to_keep = set()
    for zone in set([z for z2 in zones for z in z2.split(";")]):
        to_keep.update(construct_hull(zone,zones,g,grid_graph,adjacent_edges_crossing,grid_nodes))
    for e in grid_graph.edges():
        if not e in to_keep:
            grid_graph.delEdge(e)
    return grid_nodes



def construct_edge_regions(g,radius):
    
    return None
def insert_elements(g,grid_graph,radius,zones,grouped_data,mapping):
    
    vl = g.getLayoutProperty("viewLayout")
    vl_grid = grid_graph.getLayoutProperty("viewLayout")
    label_grid = grid_graph.getStringProperty("viewLabel")
    nodes = g.nodes()
    gene_nodes = []
    for i,z in enumerate(zones):
        gene_nodes.append({"pos":vl[nodes[i]],"nodes":[]})
    zones = {z:i for i,z in enumerate(zones)}
    
    rng = np.random.default_rng()
    for i in grouped_data.filter(["id","EnsemblID"]).itertuples():
        zone_center = vl[nodes[zones[i.id]]]
        n = grid_graph.addNode()
        zone_center = np.array((zone_center.x(),zone_center.y()))
        a= rng.uniform(0,2*np.pi)
        pos = zone_center# + np.stack((np.cos(a),np.sin(a)))*(rng.uniform()**0.5)*radius/np.sqrt(3)/2
        vl_grid[n]=tlp.Vec3f(pos[0],pos[1],0.0)
        mapping[i.Index]=n
        label_grid[n]=i.Index
        gene_nodes[zones[i.id]]["nodes"].append(n)
    return mapping,gene_nodes
def add_pathways(g,all_nodes_and_edges,zones,grouped_data):
    pathway_nodes = [n for n in all_nodes_and_edges if "is_pathways" in n["data"] and n["data"]["is_pathways"]]
    pathway_edges = {}
    for n in all_nodes_and_edges:
        if "is_pathway_edge" in n["data"] and n["data"]["is_pathway_edge"]:
            if n["data"]["target"] in pathway_edges:
                pathway_edges[n["data"]["target"]].append(n)
            else:
                pathway_edges[n["data"]["target"]]=[n]

    vl = g.getLayoutProperty("viewLayout")
    IDs = g.getStringProperty("ID")
    m = {}
    gene_to_zone = grouped_data.filter(["id","EnsemblID"]).to_dict()["id"]
    zones = {z:i for i,z in enumerate(zones)}
    nodes = g.nodes()
    fixed = g.getBooleanProperty("fixed")
    fixed.setAllNodeValue(True)
    fixed.setNodeDefaultValue(False)
    
    # vl.scale(tlp.Vec3f(10,10,1))
    avg = vl.averageEdgeLength() if g.numberOfEdges()!=0 else 20
    bl = vl.getMin()
    tr = vl.getMax()
    c = np.array(((tr.x()+bl.x())/2,(tr.y()+bl.y())/2))
    d = bl.dist(tr)
    if d==0 :
        d=20

    for i,node in enumerate(pathway_nodes):
        n = g.addNode()
        fixed[n]=False
        IDs[n]=node["data"]["id"]
        m[node["data"]["id"]]=n
        pos_neighbors = np.zeros((len(pathway_edges[node["data"]["id"]]),2))
        for j,e in enumerate(pathway_edges[node["data"]["id"]]):
            z = gene_to_zone[e["data"]["source"]]
            node_index =zones[z] 
            src = nodes[node_index]
            p = vl[src]
            pos_neighbors[j,0]=p.getX()
            pos_neighbors[j,1]=p.getY()
        p = np.mean(pos_neighbors,axis=0)
        d2 = np.linalg.norm(p-c)
        if d2<1e-2:
            vl[n]=tlp.Vec3f(d*np.cos(i*np.pi*2/len(pathway_nodes))+c[0],d*np.sin(i*np.pi*2/len(pathway_nodes))+c[1],0.0)
        else:
            p2 = c+(p-c)*d/d2
            vl[n]=tlp.Vec3f((d/100)*np.cos(i*np.pi*2/len(pathway_nodes))+p2[0],(d/100)*np.sin(i*np.pi*2/len(pathway_nodes))+p2[1],0.0)



    
    for p in pathway_edges.values():
        for e in p:
            src = e["data"]["source"]
            tgt = e["data"]["target"]
            if src in m and tgt in m:
                g.addEdge(m[src],m[tgt])
            else:
                #gene protein Edge
                assert src in gene_to_zone,f"{src} {list(gene_to_zone.keys())}"
                z = gene_to_zone[src]
                assert z in zones, f"{z} {list(zones.keys())}"
                node_index =zones[z] 
                src = nodes[node_index]
                assert tgt in m
                tgt = m[tgt]
                g.addEdge(src,tgt)
            
    return avg

def add_pathways2(g,all_nodes_and_edges,zones,grouped_data):
    pathway_nodes = [n for n in all_nodes_and_edges if "is_pathways" in n["data"] and n["data"]["is_pathways"]]
    pathway_edges = [n for n in all_nodes_and_edges if "is_pathway_edge" in n["data"] and n["data"]["is_pathway_edge"]]
    IDs = g.getStringProperty("ID")
    m = {}
    gene_to_zone = grouped_data.filter(["id","EnsemblID"]).to_dict()["id"]
    zones = {z:i for i,z in enumerate(zones)}
    nodes = g.nodes()
    fixed = g.getBooleanProperty("fixed")
    fixed.setAllNodeValue(True)
    fixed.setNodeDefaultValue(False)
    vl = g.getLayoutProperty("viewLayout")
    
    # vl.scale(tlp.Vec3f(10,10,1))
    avg = vl.averageEdgeLength() if g.numberOfEdges()!=0 else 20
    bl = vl.getMin()
    tr = vl.getMax()
    c = ((tr.x()+bl.x())/2,(tr.y()+bl.y())/2)
    d = bl.dist(tr)
    if d==0 :
        d=20

    for i,node in enumerate(pathway_nodes):
        n = g.addNode()
        fixed[n]=False
        IDs[n]=node["data"]["id"]
        m[node["data"]["id"]]=n
        vl[n]=tlp.Vec3f(d*np.cos(i*np.pi*2/len(pathway_nodes))+c[0],d*np.sin(i*np.pi*2/len(pathway_nodes))+c[1],0.0)
    
    for e in pathway_edges:
        src = e["data"]["source"]
        tgt = e["data"]["target"]
        if src in m and tgt in m:
            g.addEdge(m[src],m[tgt])
        else:
            #gene protein Edge
            z = gene_to_zone[src]
            node_index =zones[z] 
            src = nodes[node_index]
            tgt = m[tgt]
            g.addEdge(src,tgt)
            
    return avg


def spring_layout(g,edgeLength,w=0.5,num_iter=25,k=None):
    fixed = []
    pos = {}
    vl = g.getLayoutProperty("viewLayout")
    G = nx.Graph()
    fixed_prop = g.getBooleanProperty("fixed")
    nodes = g.nodes()
    m = {}
    for i,n in enumerate(nodes):
        G.add_node(i)
        m[n]=i
        if fixed_prop[n]:
            fixed.append(i)
        p = vl[n]
        pos[i]=np.array([p.x(),p.y()])
    for e in g.getEdges():
        G.add_edge(m[g.source(e)],m[g.target(e)],weight=w)
        
    pos = nx.spring_layout(G,pos=pos,fixed=fixed,iterations=num_iter,weight="weight",k=k)
    for i,n in enumerate(nodes):
        p = pos[i]
        vl[n]=tlp.Vec3f(p[0],p[1],0.0)


def add_laidout_pathways(g,grid_graph,mapping):
    fixed_prop = g.getBooleanProperty("fixed")
    vl = g.getLayoutProperty("viewLayout")
    IDs = g.getStringProperty("ID")
    vl_grid = grid_graph.getLayoutProperty("viewLayout")
    m={}
    for n in g.getNodes():
        if not fixed_prop[n]:
            n2 = grid_graph.addNode()
            vl_grid[n2]=vl[n]
            mapping[IDs[n]]=n2
            m[n]=n2
    for e in g.getEdges():
        src = g.source(e)
        tgt = g.target(e)
        if not (fixed_prop[src] or fixed_prop[tgt]):
            grid_graph.addEdge(m[src],m[tgt])
    return mapping
def move_gene_towards_protein(grid_graph,gene_node,protein_node,radius,degree):
    vl = grid_graph.getLayoutProperty("viewLayout")
    vln = vl[gene_node]
    vlp = vl[protein_node]
    pn = np.array([vln.x(),vln.y()])
    pp = np.array([vlp.x(),vlp.y()])
    d = vln.dist(vlp)
    pn = pn +((pp-pn)/d)*(radius/np.sqrt(3))/(2*degree)
    vl[gene_node]=tlp.Vec3f(pn[0],pn[1],0.0)
def randomize_pos(grid_graph,node,radius,angle):
    vl = grid_graph.getLayoutProperty("viewLayout")
    vln = vl[node]
    pn = np.array([vln.x(),vln.y()])
    pn = pn +np.array((np.cos(angle),np.sin(angle)))*(radius/np.sqrt(3))/4
    vl[node]=tlp.Vec3f(pn[0],pn[1],0.0)



def add_pathways_edges(grid_graph,mapping,all_nodes_and_edges,radius,gene_nodes):
    pathway_edges = [n for n in all_nodes_and_edges if "is_pathway_edge" in n["data"] and n["data"]["is_pathway_edge"]]
    degrees = {}
    for e in pathway_edges:
        src = e["data"]["source"]
        if src in degrees:
            degrees[src]+=1
        else:
            degrees[src]=1
    label_grid = grid_graph.getStringProperty("viewLabel")
    vl = grid_graph.getLayoutProperty("viewLayout")
    to_randomize = set()
    for i in gene_nodes:
        to_randomize.update(i["nodes"])
    for e in pathway_edges:
        tgt = mapping[e["data"]["target"]]
        src = mapping[e["data"]["source"]]
        if not grid_graph.hasEdge(src,tgt,False):
            grid_graph.addEdge(src,tgt)
            if src in to_randomize:
                move_gene_towards_protein(grid_graph,src,tgt,radius,degrees[e["data"]["source"]])
                # if label_grid[src]=="ENSG00000175899" or label_grid[src]=="ENSG00000183087":
                    # print("move",label_grid[src], vl[src])
                # to_randomize.discard(src)
    for i,n in enumerate(to_randomize):
        # if grid_graph.deg(n)==0:
            randomize_pos(grid_graph,n,radius,2*np.pi*i/len(to_randomize))
            # if label_grid[n]=="ENSG00000175899" or label_grid[n]=="ENSG00000183087":
                # print("rand", label_grid[n],vl[n])
def get_next_node(g,n,prev):
    if prev is None:
        return next(g.getInOutNodes(n))
    inout = list(g.getInOutNodes(n))
    i = inout.index(prev)
    return inout[(i+1)%len(inout)]
def get_prev_node(g,n,prev):
    inout = list(g.getInOutNodes(n))
    if prev is None:
        return inout[-1]
    i = inout.index(prev)
    return inout[(i-1)%len(inout)]
def get_hull(zones,grid_graph):
    hulls = {}
    vl = grid_graph.getLayoutProperty('viewLayout')
    for zone in set([z for z2 in zones for z in z2.split(";")]):
        sg = grid_graph.getSubGraph(f"zone_{zone}")
        tlp.PlanarityTest.planarEmbedding(sg)
        hull = []
        start = sg.nodes()[0]
        vla = sg.getStringProperty("viewLabel")
        start = get_left_most_node(sg)
        n = start
        inner = sg.getBooleanProperty("innerface")
        inner.setAllNodeValue(True)
        def visit(n):
            p = vl[n]
            hull.append((p.x(),p.y()))
            inner[n]=False
        visit(start)
        n = get_next_node(sg,start,None)
        prev =start
        while n!=start:
            visit(n)
            n2 = get_next_node(sg,n,prev)
            prev=n
            n=n2
        g2 = grid_graph.inducedSubGraph(sg.nodes())
        i=0
        while( inner.getNodesEqualTo(True).hasNext()):
            innerG = sg.inducedSubGraph(inner)
            startInner = get_left_most_node(innerG)
            path = find_path(g2,start,startInner)
            for n in path : 
                visit(n)
            n = get_prev_node(g2,startInner,path[-2])
            while not (innerG.isElement(n)):
                n = get_prev_node(g2,startInner,n)

            prev =startInner
            while n!=startInner:
                visit(n)
                n2 = get_prev_node(sg,n,prev)
                prev=n
                n=n2
            path.reverse()
            for n in path : 
                visit(n)
            sg.delSubGraph(innerG)
            i+=1
        grid_graph.delSubGraph(g2)
        hulls[zone]=hull
        


    return hulls

def find_faces(g:nx.algorithms.planarity.PlanarEmbedding):
    traversed_half_edges = dict()
    faces = []
    nodes_faces = {}
    for a in g.nodes():
        nodes_faces[a]={"cw":-1,"ccw":-1}
    for a in g.nodes():
        neighbors = g.neighbors_cw_order(a)
        for b in neighbors:
            if(a,b) not in traversed_half_edges:
                # nodes_faces[a]["cw"]=(traversed_half_edges[(a,b)])
            # else:
                face = [a,b]
                c=a
                d=b
                while d!=a:
                    e,f = g.next_face_half_edge(c,d)
                    assert e==d
                    face.append(f)
                    c=e
                    d=f
                for i in range(len(face)-1):
                    traversed_half_edges[(face[i],face[i+1])] = len(faces)
                faces.append(face)
    return traversed_half_edges,faces
def find_tlp_faces(g):
    G = nx.algorithms.planarity.PlanarEmbedding()
    mapping = {}
    G.add_nodes_from(g.nodes())
    vl = g.getLayoutProperty("viewLayout")
    vl.computeEmbedding()
    cc= tlp.ConnectedTest.computeConnectedComponents(g)
    outers = []

    left_most = [None for i in range(len(cc))]
    left = [np.inf for i in range(len(cc))]
    for i in range(len(cc)):
        for n in cc[i]:
            ref = None
            for e in g.getInOutEdges(n):
                n2= g.opposite(e,n)
                G.add_half_edge_cw(n, n2, ref)
                ref=n2
            p = vl[n]
            if p.getX()<left[i]:
                left[i] = p.getX()
                left_most[i] = n
    

        # G.add_edge(g.source(e),g.target(e))
    # traversed_half_edges,faces = find_faces(nx.check_planarity(G)[1])
    traversed_half_edges,faces = find_faces(G)
    for i in range(len(cc)):

        p1= vl[left_most[i]]
        following = None
        highest_slope = -np.inf
        for o in g.getInOutNodes(left_most[i]):
            p2 =vl[o]
            if p2.getX()==left[i]:
                slope=np.inf
            else:
                slope = (p2.getY()-p1.getY())/(p2.getX()-left[i])
            if slope > highest_slope:
                highest_slope = slope
                following = o
        outers.append(traversed_half_edges[(left_most[i],following)])
    return traversed_half_edges,faces,outers

def find_path(g,src,tgt):
    vs = g.getBooleanProperty("viewSelection")
    tlp.selectShortestPaths(g, src, tgt, tlp.OnePath, None, vs)
    path = [src]
    vs[src]=False
    while src !=tgt:
        for n in g.getInOutNodes(src):
            if vs[n]:
                path.append(n)
                vs[n]=False
                src = n
                break
    return path

def get_left_most_node(g):
    vl = g.getLayoutProperty("viewLayout")
    left_most = None
    left = np.inf
    for n in g.getNodes():
        if(vl[n].getX()<left):
            left_most=n
            left = vl[n].getX()
    return left_most


def get_surrounding_edges(grid_nodes,traversed_half_edges,faces,gene_nodes,grid_graph,outer_faces_edges,outer_faces):
    vl = grid_graph.getLayoutProperty("viewLayout")
    surrounding_edges={}
    for z in grid_nodes:
        for n in z:
            node_faces = set()
            for e in grid_graph.getInOutEdges(n):
                n2 = grid_graph.opposite(e,n)
                node_faces.add(traversed_half_edges[(n,n2)])
                node_faces.add(traversed_half_edges[(n2,n)])
            edges=set()
            is_on_outer_face = False
            for f in node_faces:
                is_on_outer_face = is_on_outer_face or f in outer_faces
                face = faces[f]
                for i in range(len(face)-1):
                    e = grid_graph.existEdge(face[i],face[i+1],False)
                    edges.add(e)
            if is_on_outer_face:
                edges.update(outer_faces_edges)
            surrounding_edges[n]=list(edges)
            # surrounding_edges[n]=grid_graph.edges()
    for i in range(len(gene_nodes)):
        j=0
        n1 = grid_nodes[i][j]
        n2 = grid_nodes[i][j+1]
        while (n1,n2) not in traversed_half_edges and j<len(grid_nodes[i]):
            j+=1
            n1 = grid_nodes[i][j]
            n2 = grid_nodes[i][(j+1)%len(grid_nodes[i])]
        f = traversed_half_edges[(n1,n2)]

        edges=set()
        face = faces[f]
        for j in range(len(face)-1):
            e = grid_graph.existEdge(face[j],face[j+1],False)
            edges.add(e)
        edges = list(edges)
        for n in gene_nodes[i]["nodes"]:
            # surrounding_edges[n] = grid_graph.edges()
            surrounding_edges[n] = edges
    return surrounding_edges


def make_surrounding_edges(grid_graph,gene_nodes,zones,grid_nodes):
    def test_tlp(graph,zone):
        edges = []
        def get_prev_node(g,n,prev):
            inout = list(g.getInOutNodes(n))
            if prev is None:
                return inout[-1]
            i = inout.index(prev)
            return inout[(i-1)%len(inout)]
        def find_next_start(g):
            zone_nodes = g["viewLabel"].getNodesEqualTo(zone)
            for start in zone_nodes:

                l = g["viewLabel"][start]
                for n in g.getOutNodes(start):
                    if l==g["viewLabel"][n]:
                        return start,n
        start,n = find_next_start(graph)

        edges.append(graph.existEdge(start,n,tlp.UNDIRECTED))
        prev = start
        
        while n!=start:
            n2 = get_prev_node(graph,n,prev)
            edges.append(graph.existEdge(n,n2,tlp.UNDIRECTED))            
            prev = n
            n=n2
        return edges
    surrounding_edges = {}
    grid_edges = set()
    grid_graph["viewLayout"].computeEmbedding()

    for i in range(len(gene_nodes)):
        zone = zones[i]

        edges = test_tlp(grid_graph,zone)
        for g in gene_nodes[i]["nodes"]:
            surrounding_edges[g]=edges
        grid_edges.update(edges)
    for z in grid_nodes:
        for n in z:
            surrounding_edges[n]=grid_edges
    return surrounding_edges,grid_edges


def euler_layout(data,all_nodes_and_edges,max_gene_size):
    # pr = cProfile.Profile()
    # pr.enable()


    # times = [time.perf_counter()]
    grouped_data = data.groupby("EnsemblID").agg(lambda x:";".join(sorted(list(set(x)))))

    g,zones,value_counts = get_planar_intersection_graph(grouped_data)
    apply_planar_layout(g)
    # tlp.saveGraph(g,"/tmp/not_grid_graph-1.tlpb.gz")
    vl = g.getLayoutProperty("viewLayout")
    # if g.numberOfNodes()>1:
    #     if g.numberOfEdges()>1:
    #         s = 1/vl.averageEdgeLength()
    #         vl.scale(tlp.Vec3f(s,s,1.0))
    #     delta = 1 #vl.averageEdgeLength() if g.numberOfEdges()>0 else 5
    #     gamma = delta/2
    #     apply_bertault(g,{n : g.edges() for n in g.nodes()} ,n_iter=30,optimalDist=delta,gamma=gamma,bends=False)

    # times.append(time.perf_counter())

    radius = compute_circle_set_radius(g)
    # print("radius", radius,vl.averageEdgeLength())
    # if np.pi*(radius**2) <np.max(list(value_counts.values()))*400:
    # n = np.max(list(value_counts.values()))
        # f = (n*400/(np.pi*(radius**2)))**0.5
    # if radius**2*np.pi/n < 1:
        # f = 1/(radius*2*np.pi/n)
        # vl.scale(tlp.Vec3f(f,f,1))
    # print("scale",f,vl.averageEdgeLength())
        # radius = compute_circle_set_radius(g)
    # print("radius", radius)
    # tlp.saveGraph(g,"/tmp/not_grid_graph0.tlpb.gz")
    # times.append(time.perf_counter())
    grid_graph = tlp.newGraph()
    grid_nodes = construct_regions(g,radius,grid_graph,zones,value_counts)
    # tlp.saveGraph(grid_graph,"/tmp/regions.tlpb.gz")

    # traversed_half_edges,faces,outer_faces = find_tlp_faces(grid_graph)
    # outer_face_edges = [grid_graph.existEdge(faces[outer_faces[j]][i],faces[outer_faces[j]][(i+1)],False) for j in range(len(outer_faces)) for i in range(len(faces[outer_faces[j]])-1)]

    vlgrid = grid_graph.getLayoutProperty("viewLayout")
    # print("a",vlgrid.averageEdgeLength(),tlp.computeBoundingBox(grid_graph))
    # apply_bertault(grid_graph)
    # tlp.saveGraph(grid_graph,"/tmp/grid_graph1.tlpb.gz")
    # times.append(time.perf_counter())

    avg = add_pathways(g,all_nodes_and_edges,zones,grouped_data)
    # print("b",vl.averageEdgeLength())
    # tlp.saveGraph(g,"/tmp/not_grid_graph1.tlpb.gz")
    # times.append(time.perf_counter())
    # spring_layout(g,avg,w=1,k=avg)
    # print("c",vl.averageEdgeLength())
    # tlp.saveGraph(g,"/tmp/not_grid_graph2.tlpb.gz")
    fixed = grid_graph.getBooleanProperty("fixed")
    fixed.setValueToGraphNodes(True,grid_graph)
    fixed.setNodeDefaultValue(False)
    grid_size=grid_graph.numberOfNodes()
    # times.append(time.perf_counter())
    
    mapping = add_laidout_pathways(g,grid_graph,{})
    pathway_nodes = set(mapping.values())
    # print("d",vlgrid.averageEdgeLength())
    # tlp.saveGraph(grid_graph,"/tmp/grid_graph2.tlpb.gz")
    fixed.setNodeDefaultValue(True)
    # print("dlghdlfgdklfgn" ,grid_graph.numberOfNodes())
    # times.append(time.perf_counter())
    mapping,gene_nodes = insert_elements(g,grid_graph,radius/10,zones,grouped_data,mapping)
    # times.append(time.perf_counter())
    surrounding_edges,grid_edges =  make_surrounding_edges(grid_graph,gene_nodes,zones,grid_nodes)
    # times.append(time.perf_counter())

    # print("dlghdlfgdklfgn" ,grid_graph.numberOfNodes())
    # surrounding_edges = get_surrounding_edges(grid_nodes,traversed_half_edges,faces,gene_nodes,grid_graph,outer_face_edges,outer_faces)
    # print("e",vlgrid.averageEdgeLength())
    # tlp.saveGraph(grid_graph,"/tmp/grid_graph3.tlpb.gz")

    # print("f",vlgrid.averageEdgeLength())
    # test = grid_graph.inducedSubGraph(only_grid.nodes()+[n for z in gene_nodes for n in z["nodes"]])
    # vl_test = test.getLayoutProperty("viewLayout")
    vla = grid_graph.getStringProperty("viewLabel")
    # vla_test = test.getStringProperty("viewLabel")
    # to_remove = []
    # # TODO add edges
    # # for i,z in enumerate(gene_nodes):
    # #     c = test.addNode()
    # #     to_remove.append(c)
    # #     vl_test[c]=z["pos"]
    # #     vla_test[c]=f"center_zone_{i}"
    #     # for n in z["nodes"]:
    #         # test.addEdge(n,c)
    # tlp.saveGraph(test,"/tmp/test0.tlpb.gz")
    # # apply_bertault(test)
    
    # apply_planar_layout(test,(vl_test.getMax().getX()-vl_test.getMin().getX())/2)
    # tlp.saveGraph(test,"/tmp/test.tlpb.gz")
    # test.delNodes(to_remove)
    only_grid = grid_graph.inducedSubGraph(grid_graph.nodes())
    add_pathways_edges(grid_graph,mapping,all_nodes_and_edges,radius/5,gene_nodes)
    avg_only_grid = vlgrid.averageEdgeLength(only_grid)
    # print("onlygrid_edge_length",avg_only_grid)
    # if(only_grid.numberOfNodes()>100):
    #     print(vlgrid[only_grid.nodes()[98]],vlgrid[only_grid.nodes()[100]])
    #     print(vla[only_grid.nodes()[98]],vla[only_grid.nodes()[100]])
    vlgrid.center(tlp.computeBoundingRadius(only_grid)[0],grid_graph)
    # outer_face = grid_graph.getBooleanProperty("outerface")
    # outer_face.setAllNodeValue(False)
    # outer_face.setAllEdgeValue(False)
    # for e in outer_face_edges:
    #     outer_face[e]=True
    #     outer_face[grid_graph.source(e)]=True
    #     outer_face[grid_graph.target(e)]=True
    # tlp.saveGraph(grid_graph,"/tmp/grid_graph6.tlpb.gz")
    # times.append(time.perf_counter())
    apply_bertault(only_grid,surrounding_edges,n_iter=25,gamma=avg_only_grid/2,optimalDist=avg_only_grid,bends=False,flexible_edges=False,fixed_nodes=pathway_nodes,debug=True)
    # times.append(time.perf_counter())
    # tlp.saveGraph(grid_graph,"/tmp/grid_graph7.tlpb.gz")
    # apply_bertault(grid_graph,surrounding_edges,n_iter=5,gamma=avg_only_grid/8,optimalDist=avg_only_grid,bends=False,flexible_edges=False,fixed_nodes=set(grid_graph.nodes()).difference(mapping.values()))
    # tlp.saveGraph(grid_graph,"/tmp/grid_graph4.tlpb.gz")
    for n in pathway_nodes:
        # surrounding_edges[n]=outer_face_edges
        surrounding_edges[n]=grid_edges
        # edges = set()
        # for n2 in grid_graph.getInOutNodes(n):
        #     edges.update(surrounding_edges[n2])
        # surrounding_edges[n]=edges

    # tlp.saveGraph(grid_graph,"/tmp/grid_graphtest.tlpb.gz")
    # print(avg_only_grid,vlgrid.averageEdgeLength(only_grid))
    vlgrid.center(grid_graph)
    apply_bertault(grid_graph,surrounding_edges,n_iter=20,fixed_nodes = set([n for n in grid_graph.nodes() if n not in pathway_nodes]),gamma=avg_only_grid*2,optimalDist=avg_only_grid*2,bends=False,flexible_edges=False,debug=False)
    # times.append(time.perf_counter())
    # apply_bertault(grid_graph)
    # tlp.saveGraph(grid_graph,"/tmp/grid_graph6.tlpb.gz")
    avg = vlgrid.averageEdgeLength()
    
    vl = grid_graph.getLayoutProperty("viewLayout")
    not_pathway_nodes = list(set(grid_graph.nodes()).difference(pathway_nodes))
    minDist = np.inf
    for i,u in enumerate(not_pathway_nodes):
            for j in range(i+1,len(not_pathway_nodes)):
                    v = not_pathway_nodes[j]
                    d = vl[u].dist(vl[v])
                    if d<minDist:
                        minDist = d
                        # print(vla[u],vla[v],u,v,minDist)
                    
    # print("minDist<max_gene_size*10:",minDist,max_gene_size,"*10",vl.averageEdgeLength())

    if minDist<np.sqrt(max_gene_size)*20:
        vl.scale(tlp.Vec3f(np.sqrt(max_gene_size)*20/minDist,np.sqrt(max_gene_size)*20/minDist,1.0))
    elif minDist>np.sqrt(max_gene_size)*20*3:
        vl.scale(tlp.Vec3f(np.sqrt(max_gene_size)*20*3/minDist,np.sqrt(max_gene_size)*20*3/minDist,1.0))
    # vl.scale(tlp.Vec3f(100.0,100.0,0.0))
    # r = tlp.computeBoundingRadius(grid_graph)
    # print(r[1].dist(r[0]))
    pos = {}
    for n in all_nodes_and_edges:
        if "source" not in n["data"]:
            p = vl[mapping[n["data"]["id"]]]
            pos[n["data"]["id"]]=(p.x(),p.y())
            n["position"]={'x':p.x(),'y':p.y()}
        pos = {}
    for n in all_nodes_and_edges:
        if "source" not in n["data"]:
            p = vl[mapping[n["data"]["id"]]]
            pos[n["data"]["id"]]=(p.x(),p.y())
            n["position"]={'x':p.x(),'y':p.y()}
    # times.append(time.perf_counter())

    hull  = get_hull(zones,grid_graph)
    # times.append(time.perf_counter())
    # times = np.array(times)
    # np.set_printoptions(linewidth=750)
    # print(times[1:]-times[:-1])
    # print((times[1:]-times[:-1])/(times[-1]-times[0]))
    # pr.disable()
    # s = io.StringIO()
    # sortby = pstats.SortKey.CUMULATIVE
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats(10)
    # print(s.getvalue())
    return pos,hull



if __name__ == "__main__":
    import detail_graph
    import data_manager
    import sys
    import tqdm
    dm = data_manager.DataManager(sys.argv[1],sys.argv[2],sys.argv[3])

    with open("/tmp/test_signs2/part-00000") as f:
        for l in tqdm.tqdm(f.readlines()):
            l = [ i for i in l[:-1].split(";") if len(i)>0]
            # if len(l)<=12:
            updated_elements = []
            edges = []
            nodes = []
            stylesheet_detail = []

            # print(updated_elements)
            # print(existing_elements)
            # updated_elements.clear()
            intersections,items = dm.get_genes_intersections(id_filter=l,disease_filter=["KIRC","HNSC","LUAD","LIHC"],gene_filter=[],selected_filter="Merge",comparisons_filter=["NormalvsStageI","NormalvsStageII","NormalvsStageIII","NormalvsStageIV"])
            data = dm.get_genes_intersection_data(id_filter=l,disease_filter=["KIRC","HNSC","LUAD","LIHC"],gene_filter=[],selected_filter="Merge",comparisons_filter=["NormalvsStageI","NormalvsStageII","NormalvsStageIII","NormalvsStageIV"])


            sizes = items.set_index("gene")
            updated = set()
            # update_cur_elements(existing_elements,updated,stylesheet_detail,sizes,updated_elements)
            symbols = dm.get_symbol([g.gene for g in items.itertuples() if g.gene not in updated])
            all_nodes = [{'data':{'id':g.gene,"label":"","weight":(g.size),"Signatures":g.id,"tooltip_content":symbols[g.gene]}} for g in items.itertuples()]
            nodes = [{'data':{'id':g.gene,"label":"","weight":(g.size),"Signatures":g.id,"tooltip_content":symbols[g.gene]}} for g in items.itertuples() if g.gene not in updated]

            all_nodes_and_edges,stylesheet_detail = detail_graph.add_path_ways(all_nodes,stylesheet_detail,updated_elements,True,dm)

            grouped_data = data.groupby("EnsemblID").agg(lambda x:";".join(sorted(list(set(x)))))

            pos,hull = euler_layout(data,all_nodes_and_edges,data.filter(["EnsemblID","id"]).groupby("EnsemblID").agg(lambda a: len(a))["id"].max())
            hull_polygons = {z:shapely.polygons(h) for z,h in hull.items()}
            for p in hull_polygons.values():
                shapely.prepare(p)
            for i in grouped_data.filter(["id","EnsemblID"]).itertuples():
                p = pos[i.Index]
                for z in i.id.split(";"):
                    assert shapely.contains_xy(hull_polygons[z],*p),(z,i,pos)
            test_pos = np.array(list(pos.values()))
            bb = (np.min(test_pos,axis=0),np.max(test_pos,axis=0))
            diag = ((bb[1][0]-bb[0][0])**2+(bb[1][1]-bb[0][1])**2)**0.5
            assert diag<10000,(diag,l)


