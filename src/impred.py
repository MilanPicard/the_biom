import time
import numpy as np
import tqdm
from tulip import tlp
from scipy import spatial as sspatial
# from py4j.java_gateway import JavaGateway
import pandas as pd
class Impred:
    def __init__(self,g,delta,gamma,nb_iter,zones,surrounding_edges,fixed_nodes) -> None:
        self.g = g

        self.delta = delta
        self.gamma = gamma
        self.nb_iter = nb_iter
        self.vl = g.getLayoutProperty("viewLayout")

        self.center = tlp.computeBoundingBox(g).center()
        # self.vl.center()
        # avg = self.vl.averageEdgeLength()
        # self.scale = avg/self.delta
        self.scale = 1
        # self.vl.scale(tlp.Vec3f(1/self.scale,1/self.scale,1.0))
        pos = []
        for n in g.nodes():
            p = self.vl[n]
            pos.append((p.getX(),p.getY()))
        self.pos = np.array(pos)
        self.kdt = sspatial.KDTree(self.pos)
        self.zones = zones
        adj = []
        for e in g.getEdges():
            adj.append((g.nodePos(g.source(e)),g.nodePos(g.target(e))))
        if len(adj)==0:
            self.adj = np.empty((0,2),dtype=int)
        else:
            self.adj = np.array(adj)
        # self.gateway = JavaGateway()
        self.surrounding_edges =surrounding_edges
        self.fixed_nodes = fixed_nodes
        self.movable_nodes = list(set(g.nodes()).difference(fixed_nodes))
        self.movable_index = [g.nodePos(n) for n in self.movable_nodes]
        self.is_movable = np.zeros(len(self.pos),dtype=bool)
        self.is_movable[self.movable_index]=True
        point_index, src_index, tgt_index = self.get_edge_rep_indices()
        self.edge_rep_point_index = point_index
        self.edge_rep_src_index=src_index
        self.edge_rep_tgt_index= tgt_index


    def layout(self,debug=False):
        self.compute_faces()
        global_maxes = np.linspace(3*self.delta,0,num=self.nb_iter)
        vl=self.g.getLayoutProperty("viewLayout")
        for j,global_max in enumerate(global_maxes):
            self.iter(global_max,debug)
        for i,n in enumerate(self.g.nodes()):
            vl[n]=tlp.Vec3f(self.pos[i,0],self.pos[i,1],0.0)
        
        vl.scale(tlp.Vec3f(self.scale,self.scale,1.0))
        vl.center(self.center)
        return self.pos + np.array([[self.center.getX(),self.center.getY()]])
    def iter(self,global_max,debug=False):
        f = self.compute_forces()
        max_move=self.compute_max_move(global_max,f)
        self.move(f,max_move)
    def compute_forces(self,):
        fg = self.gravity_force()
        fn = self.node_node_rep()
        fa = self.edge_attr()
        fe = self.edge_rep()
        r = fn+fa+fe+fg
        return r
    def gravity_force(self):
        f = np.zeros_like(self.pos)
        diff = -self.pos[self.movable_index]
        f[self.movable_index] = diff*(np.sum(diff**2,-1,keepdims=True))
        return f

    def node_node_rep(self):
        pairs = self.kdt.query_pairs(r=3*self.delta,output_type="ndarray")
        first_movable = self.is_movable[pairs[:,0]]
        second_movable = self.is_movable[pairs[:,1]]
        one_movable = np.logical_or(first_movable,second_movable)
        forces = np.zeros_like(self.pos)
        if np.any(one_movable):
            pairs = pairs[one_movable]
            first_movable = first_movable[one_movable]
            second_movable = second_movable[one_movable]
            diff = self.pos[pairs[:,0]]-self.pos[pairs[:,1]]
            dist = np.linalg.norm(diff,axis=1)
            mag = ((self.delta)/dist)**2/dist

            repulsions = np.expand_dims(mag,-1) *diff
            forces[:,0]+=np.bincount(pairs[first_movable,0],repulsions[first_movable,0],minlength=self.pos.shape[0])
            forces[:,1]+=np.bincount(pairs[first_movable,0],repulsions[first_movable,1],minlength=self.pos.shape[0])
            forces[:,0]+=np.bincount(pairs[second_movable,1],-repulsions[second_movable,0],minlength=self.pos.shape[0])
            forces[:,1]+=np.bincount(pairs[second_movable,1],-repulsions[second_movable,1],minlength=self.pos.shape[0])

        return forces
    def edge_attr(self):
        first_movable = self.is_movable[self.adj[:,0]]
        second_movable = self.is_movable[self.adj[:,1]]
        one_movable = np.logical_or(first_movable,second_movable)
        adj = self.adj[one_movable]
        first_movable = first_movable[one_movable]
        second_movable = second_movable[one_movable]

        diff = self.pos[adj[:,0]]-self.pos[adj[:,1]]
        dist = np.linalg.norm(diff,axis=1)
        mag = ((dist/self.gamma)/self.gamma)#**0.4
        attractions = np.expand_dims(mag,-1) *diff
        forces = np.zeros_like(self.pos)
        forces[:,0]+=np.bincount(adj[second_movable,1],attractions[second_movable,0],minlength=self.pos.shape[0])

        forces[:,1]+=np.bincount(adj[second_movable,1],attractions[second_movable,1],minlength=self.pos.shape[0])
        forces[:,0]+=np.bincount(adj[first_movable,0],-attractions[first_movable,0],minlength=self.pos.shape[0])
        forces[:,1]+=np.bincount(adj[first_movable,0],-attractions[first_movable,1],minlength=self.pos.shape[0])
        return forces
    def edge_rep(self):
        # indices = self.kdt.query_ball_point(self.pos,self.gamma)
        # nodes = self.g.nodes()
        point_index, src_index, tgt_index = self.edge_rep_point_index,self.edge_rep_src_index,self.edge_rep_tgt_index

        forces = np.zeros_like(self.pos)
        pos_src = []
        pos_tgt = []
        pos_point = []
        if(len(src_index)!=0):
            pos_src = self.pos[src_index]
            pos_tgt = self.pos[tgt_index]
            pos_point = self.pos[point_index]
            proj = self.get_projections(pos_src,pos_tgt,pos_point)
            proj_dist = np.linalg.norm(proj,axis=1)
            close_enough = proj_dist<self.gamma*3
            pos_src3 = pos_src[close_enough]
            pos_tgt3 = pos_tgt[close_enough]
            pos_point3 = pos_point[close_enough]
            point_index3 = point_index[close_enough]
            tgt_index3 = tgt_index[close_enough]
            src_index3 = src_index[close_enough]
            proj3 = proj[close_enough]
            # print(len(close_enough),len(proj))
            on_edge = np.logical_and(np.logical_or(np.abs(pos_src3[:,0]-pos_tgt3[:,0])<1e-4,np.logical_xor(pos_src3[:,0]<=proj3[:,0],pos_tgt3[:,0]<=proj3[:,0])),
                                    np.logical_or(np.abs(pos_src3[:,1]-pos_tgt3[:,1])<1e-4,np.logical_xor(pos_src3[:,1]<=proj3[:,1],pos_tgt3[:,1]<=proj3[:,1])))
            if np.any(on_edge):
                pos_point2 = pos_point3[on_edge]
                pos_src2 = pos_src3[on_edge]
                pos_tgt2 = pos_tgt3[on_edge]
                point_index2 = point_index3[on_edge]
                tgt_index2 = tgt_index3[on_edge]
                src_index2 = src_index3[on_edge]
                proj2 = proj3[on_edge]
                diff = pos_point2-proj2
                dist = np.linalg.norm(diff,axis=1)
                # cond = dist<self.gamma
                # proj = proj[cond]
                # pos_point = pos_point[cond]
                # point_index = point_index[cond]
                mag = (self.gamma/dist)**2/dist
                repulsion = np.expand_dims(mag,-1)*(diff)
                # for i in range(len(pos_point2)):
                    # print("inedge",pos_point2[i],pos_src2[i],pos_tgt2[i],proj2[i],repulsion[i],dist[i],self.gamma,"my",sep=";")
                forces[:,0]+=np.bincount(point_index2,repulsion[:,0],minlength=self.pos.shape[0])
                forces[:,1]+=np.bincount(point_index2,repulsion[:,1],minlength=self.pos.shape[0])
                balance = np.linalg.norm(proj2-pos_tgt2,axis=1)/np.linalg.norm(pos_src2-pos_tgt2,axis=1)
                src_movable = self.is_movable[src_index2]
                tgt_movable = self.is_movable[tgt_index2]
                if np.any(src_movable):
                    forces[:,0]-=np.bincount(src_index2[src_movable],repulsion[src_movable,0]*balance,minlength=self.pos.shape[0])
                    forces[:,1]-=np.bincount(src_index2[src_movable],repulsion[src_movable,1]*balance,minlength=self.pos.shape[0])
                if np.any(tgt_movable):
                    forces[:,0]-=np.bincount(tgt_index2[tgt_movable],repulsion[tgt_movable,0]*(1-balance),minlength=self.pos.shape[0])
                    forces[:,1]-=np.bincount(tgt_index2[tgt_movable],repulsion[tgt_movable,1]*(1-balance),minlength=self.pos.shape[0])
            if not np.all(on_edge):
                pos_point2 = pos_point3[np.logical_not(on_edge)]
                pos_src2 = pos_src3[np.logical_not(on_edge)]
                pos_tgt2 = pos_tgt3[np.logical_not(on_edge)]
                point_index2 = point_index3[np.logical_not(on_edge)]
                tgt_index2 = tgt_index3[np.logical_not(on_edge)]
                src_index2 = src_index3[np.logical_not(on_edge)]
                proj2 = proj3[np.logical_not(on_edge)]
                sdiff = pos_point2-pos_src2
                tdiff = pos_point2-pos_tgt2
                sdist = np.linalg.norm(sdiff,axis=1)
                tdist = np.linalg.norm(tdiff,axis=1)
                # cond = dist<self.gamma
                # proj = proj[cond]
                # pos_point = pos_point[cond]
                # point_index = point_index[cond]
                mag = (self.gamma/sdist)**2/sdist
                repulsion = np.expand_dims(mag,-1)*(sdiff)
                # for i in range(len(pos_point2)):
                #     print("not inedge",pos_point2[i],pos_src2[i],pos_tgt2[i],proj2[i],repulsion[i],"my",sep=";")
                forces[:,0]+=np.bincount(point_index2,repulsion[:,0],minlength=self.pos.shape[0])
                forces[:,1]+=np.bincount(point_index2,repulsion[:,1],minlength=self.pos.shape[0])
                # balance = np.linalg.norm(proj2-pos_tgt2,axis=1)/np.linalg.norm(pos_src2-pos_tgt2,axis=1)
                src_movable = self.is_movable[src_index2]
                tgt_movable = self.is_movable[tgt_index2]
                if np.any(src_movable):
                    forces[:,0]-=np.bincount(src_index2[src_movable],repulsion[src_movable,0],minlength=self.pos.shape[0])
                    forces[:,1]-=np.bincount(src_index2[src_movable],repulsion[src_movable,1],minlength=self.pos.shape[0])
                mag = (self.gamma/tdist)**2/tdist
                repulsion = np.expand_dims(mag,-1)*(tdiff)
                forces[:,0]+=np.bincount(point_index2,repulsion[:,0],minlength=self.pos.shape[0])
                forces[:,1]+=np.bincount(point_index2,repulsion[:,1],minlength=self.pos.shape[0])
                if np.any(tgt_movable):
                    forces[:,0]-=np.bincount(tgt_index2[tgt_movable],repulsion[tgt_movable,0],minlength=self.pos.shape[0])
                    forces[:,1]-=np.bincount(tgt_index2[tgt_movable],repulsion[tgt_movable,1],minlength=self.pos.shape[0])
        return forces

    def get_edge_rep_indices(self):
        point_index = []
        src_index = []
        tgt_index = []
        for i,n in enumerate(self.movable_nodes):
            srcs = []
            tgts = []
            if n in self.surrounding_edges:
                for e in self.surrounding_edges[n]:
                    s = self.g.source(e)
                    t = self.g.target(e)
                    if not (n == s or n == t):
                        srcs.append(self.g.nodePos(s))
                        tgts.append(self.g.nodePos(t))
                if len(srcs)!=0:
                    # pos_src.append(self.pos[srcs])
                    # pos_tgt.append(self.pos[tgts])
                    j = self.g.nodePos(n)
                    # pos_point.append(np.repeat(np.expand_dims(self.pos[j],0),len(srcs),0))
                    point_index.append(np.full(len(srcs),j))
                    src_index.append(srcs)
                    tgt_index.append(tgts)
        if len(src_index)>0:
            src_index = np.concatenate(src_index)
            tgt_index = np.concatenate(tgt_index)
            point_index = np.concatenate(point_index)
        return point_index,src_index,tgt_index

        
    def get_projections(self,pos_src,pos_tgt,pos_point,debug=False):
        v = pos_tgt-pos_src
        d = np.linalg.norm(v,axis=-1)
        unit = v/np.expand_dims(d,-1)
        dot = np.sum( (pos_point-pos_src)* unit,-1,keepdims=True)
        # dot2 = np.stack([np.dot(pos_point[i]-pos_src[i], unit[i]) for i in range(len(pos_src))],0)
        # assert(np.all(dot==dot2)), np.max(np.abs(dot-dot2))
        r= unit * dot+pos_src
        # if debug:
        #     for i in range(len(r)):
        #         tgt = self.gateway.jvm.ocotillo.geometry.Coordinates(pos_tgt[i,0],pos_tgt[i,1],0.0)
        #         src = self.gateway.jvm.ocotillo.geometry.Coordinates(pos_src[i,0],pos_src[i,1],0.0)
        #         p = self.gateway.jvm.ocotillo.geometry.Coordinates(pos_point[i,0],pos_point[i,1],0.0)
        #         v2 = tgt.minus(src)
        #         assert np.linalg.norm(v[i]-np.array([v2.x(),v2.y()]))<1e-6,(np.linalg.norm(v[i]-np.array([v2.x(),v2.y()])),v[i],v2.toString())
        #         l = self.gateway.jvm.ocotillo.geometry.Geom2D.magnitude(v2)
        #         assert np.abs(l-d[i])<1e-6,(np.abs(l-d[i]),d[i],l)
        #         unit2 = v2.divide(l)
        #         assert np.linalg.norm(unit[i]-np.array([unit2.x(),unit2.y()]))<1e-6,(np.linalg.norm(unit[i]-np.array([unit2.x(),unit2.y()])),unit[i],unit2.toString())
        #         dot2 = self.gateway.jvm.ocotillo.geometry.Geom2D.dotProduct(p.minus(src),unit2)
        #         assert np.abs(dot2-dot[i])<1e-6,(np.abs(dot2-dot[i]),dot[i],dot2)
        #         r2 = unit2.times(dot2).plus(src)
        #         assert np.linalg.norm(r[i]-np.array([r2.x(),r2.y()]))<1e-6,(np.linalg.norm(r[i]-np.array([r2.x(),r2.y()])),r[i],r2.toString())
        #         print(r,r2.toString(),self.gateway.jvm.ocotillo.geometry.Geom2D.pointOnLineProjection(p,src,tgt).toString())
        return r

        # edge_slopes = (pos_src[:,1]-pos_tgt[:,1])/(pos_src[:,0]-pos_tgt[:,0])
        # perp_slopes = -1/edge_slopes
        # a = pos_src[:,1]-edge_slopes*pos_src[:,0]
        # b = pos_point[:,1]-perp_slopes*pos_point[:,0]
        # x = (b-1)/(edge_slopes-perp_slopes)
        # y = edge_slopes*x+b
        # return np.stack([x,y],axis=-1)
    def  test_angle(self,y,x,debug=False):
        a = np.arctan2(y , x)
        if(a>np.pi):
            a -= np.pi*2
        return a
        # if (x > 0.0 and y >= 0.0):
        #     if(debug):
        #         print("a")
        #     angle = np.arctan2(y , x)
        # elif (x > 0.0 and y < 0.0):
        #     if(debug):
        #         print("b")
        #     angle = np.arctan2(y , x) + np.pi*2;
        # elif (x < 0.0):
        #     if(debug):
        #         print("c",y,x,np.arctan2(y,x),np.arctan2(y,x)*180/np.pi)
        #     angle = np.arctan2(y , x) +np.pi;
        # elif (np.abs(x)<1e-3 and y > 0.0):
        #     if(debug):
        #         print("d")
        #     angle = np.pi/2
        # elif (np.abs(x) and  y < 0.0):
        #     if(debug):
        #         print("e")
        #     angle = np.pi*3/2
        # return angle
    def set_max_move(self,f,max_move,i,angle,distance,debug=False):
        cond1 = np.linalg.norm(f[i],axis=1)>distance
        i2 = i[cond1]
        do_something = np.any(cond1)
        if do_something:
            fa = np.arctan2(f[i2,1],f[i2,0])
            fa[fa>np.pi]-=np.pi*2
            fa = fa%(np.pi*2)
            angle_diff = (fa-angle[cond1])%(np.pi*2)
            angle_diff = np.abs(np.where(angle_diff>np.pi,angle_diff-np.pi*2,angle_diff))
            cond = np.logical_or(angle_diff<np.pi/2,angle_diff>np.pi*3/2)
            do_something = np.any(cond)

            if do_something:
                d = distance[cond1][cond]
                i3 = i2[cond]
                order = np.argsort(i3)
                i3 = i3[order]
                v = (d/np.cos(angle_diff[cond]))[order]
                uniques = np.unique(i3,return_index=True)
                v=np.array([np.min(j) for j in np.split(v,uniques[1][1:])])
                max_move[uniques[0]]=np.minimum(max_move[uniques[0]],v)

    

    def compute_max_move(self,global_max,f,debug=False):
        tz1 = time.perf_counter()
        def test(a,b,c,debug=False):
            # if debug:
            #     print(f"""
            #         {(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])}\n
            #         {(b[0]-a[0])*(c[1]-a[1])}-{(b[1]-a[1])*(c[0]-a[0])}\n
            #         {(b[0]-a[0])}*{(c[1]-a[1])}-{(b[1]-a[1])}*{(c[0]-a[0])}\n
            #     """)
            return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])>0
        
        
        max_move = np.full((self.g.numberOfNodes()),global_max)
        srcs= []
        tgts= []
        points = []
        points, srcs, tgts = self.edge_rep_point_index,self.edge_rep_src_index,self.edge_rep_tgt_index
        if len(srcs)!=0:
            pos_srcs = self.pos[srcs]
            pos_tgts = self.pos[tgts]
            # pos = np.expand_dims(self.pos[i],0).repeat(len(pos_srcs),0)
            pos = self.pos[points]
            projs = self.get_projections(pos_srcs,pos_tgts,pos)
            dist_p = np.linalg.norm(projs-pos,axis=1)
            not_too_far = dist_p<=global_max*2
            # if i==249 and not np.any(np.logical_and(
            #     np.logical_xor(projs[:,0][not_too_far]<pos_srcs[:,0][not_too_far],projs[:,0][not_too_far]<pos_tgts[not_too_far][:,0]),
            #     np.logical_xor(projs[:,1][not_too_far]<pos_srcs[:,1][not_too_far],projs[:,1][not_too_far]<pos_tgts[not_too_far][:,1]))):
            #     print("no _on _edge",not_too_far,projs,pos,pos_srcs,pos_tgts,np.logical_and(
            #     np.logical_xor(projs[:,0][not_too_far]<pos_srcs[:,0][not_too_far],projs[:,0][not_too_far]<pos_tgts[not_too_far][:,0]),
            #     np.logical_xor(projs[:,1][not_too_far]<pos_srcs[:,1][not_too_far],projs[:,1][not_too_far]<pos_tgts[not_too_far][:,1])))

            srcs = srcs[not_too_far]
            tgts = tgts[not_too_far]
            points = points[not_too_far]
            pos_srcs = pos_srcs[not_too_far]
            pos_tgts = pos_tgts[not_too_far]
            pos = pos[not_too_far]
            projs = projs[not_too_far]
            dist_p = dist_p[not_too_far]
            dist_s = np.linalg.norm(pos_srcs-pos,axis=1)
            dist_t = np.linalg.norm(pos_tgts-pos,axis=1)
            on_edge = np.logical_and(
                np.logical_or(np.abs(pos_tgts[:,0]-pos_srcs[:,0])<1e-4,np.logical_xor(projs[:,0]<pos_srcs[:,0],projs[:,0]<pos_tgts[:,0])),
                np.logical_or(np.abs(pos_tgts[:,1]-pos_srcs[:,1])<1e-4,np.logical_xor(projs[:,1]<pos_srcs[:,1],projs[:,1]<pos_tgts[:,1]))
            )

            if(np.any(on_edge)):
                # point_coord =self.gateway.jvm.ocotillo.geometry.Coordinates(pos[0,0],pos[0,1],0.0)
                on_edge_projs = projs[on_edge]
                on_edge_pos = pos[on_edge]
                angle = np.arctan2(on_edge_projs[:,1]-on_edge_pos[:,1],on_edge_projs[:,0]-on_edge_pos[:,0])%(2*np.pi)
                on_edge_dist = dist_p[on_edge]


                self.set_max_move(f,max_move,points[on_edge],angle,np.clip(on_edge_dist/2.1-self.gamma/2,0,None),debug=False)
                self.set_max_move(f,max_move,srcs[on_edge],angle+np.pi,np.clip(on_edge_dist/2.1-self.gamma/2,0,None),debug=False)
                self.set_max_move(f,max_move,tgts[on_edge],angle+np.pi,np.clip(on_edge_dist/2.1-self.gamma/2,0,None),debug=False)
                # if i==21 and len(self.pos)>21:
                    # print("debug HSNC NvsIII a",self.pos[21],f[21],max_move[21],global_max)

                # times2[6]+=ty3-ty2
                # times2[3]+=tx5-tx4
            if(not np.all(on_edge)):
                
                not_on_edge = np.logical_not(on_edge)
                pos= pos[not_on_edge]
                close = np.where(dist_s[not_on_edge]<dist_t[not_on_edge],srcs[not_on_edge],tgts[not_on_edge])
                far = np.where(dist_s[not_on_edge]>=dist_t[not_on_edge],srcs[not_on_edge],tgts[not_on_edge])
                # assert np.all(far!=close)
                close_pos = self.pos[close]
                far_pos = self.pos[far]

                angle = np.arctan2(close_pos[:,1]-pos[:,1],close_pos[:,0]-pos[:,0])%(2*np.pi)
                close_dist=np.linalg.norm(close_pos-pos,axis=-1)
                middle = (close_pos+pos)/2
                other = middle+np.stack([np.cos(angle+np.pi/2),np.sin(angle+np.pi/2)],-1)
                projs2 = self.get_projections(middle,other,far_pos)                
                # for j in range(len(pos)):
                # #     cp = self.gateway.jvm.ocotillo.geometry.Coordinates(close_pos[j,0],close_pos[j,1])
                # #     fp = self.gateway.jvm.ocotillo.geometry.Coordinates(far_pos[j,0],far_pos[j,1])
                # #     pp = self.gateway.jvm.ocotillo.geometry.Coordinates(pos[j,0],pos[j,1])
                # #     nodeAngle = self.gateway.jvm.ocotillo.geometry.Geom2D.angle(cp.minus(pp))
                # #     assert np.abs(nodeAngle-angle[j])<1e-3
                # #     edgeAngle = self.gateway.jvm.ocotillo.geometry.Geom2D.posNormalizeRadiansAngle(nodeAngle + self.gateway.jvm.Math.PI)
                # #     axisAngle = nodeAngle + 1.5707963267948966
                # #     assert np.abs(nodeAngle-angle[j])<1e-3
                # #     assert np.abs(edgeAngle-(angle[j]+np.pi)%(np.pi*2))<1e-3,(edgeAngle,(angle[j]+np.pi))
                # #     axisPointA = pp.plus(cp.minus(pp).divide(2.0))
                # #     axisPointB = self.gateway.jvm.ocotillo.geometry.Geom2D.unitVector(axisAngle).plusIP(axisPointA)
                # #     farExtrCollDist = self.gateway.jvm.ocotillo.geometry.Geom2D.pointToLineDistance(fp, axisPointA, axisPointB)
                # #     farExtrCollDist2 = np.linalg.norm(projs2-far_pos,axis=1)[j]
                # #     assert np.abs(farExtrCollDist-farExtrCollDist2)<1e-3

                #     middle2 = (pos[j]+close_pos[j])/2
                #     middle_l = middle2 + np.stack((np.cos(angle[j]+np.pi/2),np.sin(angle[j]+np.pi/2)))
                #     # t1 = test(pos_srcs[on_edge][j],pos_tgts[on_edge][j],[j])
                #     tp = test(middle2,middle_l,pos[j])
                #     ts = test(middle2,middle_l,close_pos[j])
                #     tt = test(middle2,middle_l,far_pos[j])
                #     self.set_max_move(f,max_move,np.full((1,),i),angle[j:j+1],close_dist[j:j+1]/2.1)
                #     self.set_max_move(f,max_move,close[j:j+1],angle[j:j+1]+np.pi,close_dist[j:j+1]/2.1)
                #     self.set_max_move(f,max_move,far[j:j+1],angle[j:j+1]+np.pi,np.linalg.norm(projs2[j:j+1]-far_pos[j:j+1],axis=1)/1.05)
                #     mag = np.linalg.norm(f[i],keepdims=True)
                #     if mag[0]>max_move[i]:
                #         f1 = f[i]*max_move[i]/mag
                #     else:
                #         f1 = f[i]
                #     c = self.pos[i]+f1
                #     mag = np.linalg.norm(f[close[j]],keepdims=True)
                #     if mag[0]>max_move[close[j]]:
                #         f2 = f[close[j]]*max_move[close[j]]/mag
                #     else:
                #         f2 = f[close[j]]
                #     a = close_pos[j]+f2
                    
                #     mag = np.linalg.norm(f[far[j]],keepdims=True)
                #     if mag[0]>max_move[far[j]]:
                #         f3 = f[far[j]]*max_move[far[j]]/mag
                #     else:
                #         f3 = f[far[j]]
                #     b = far_pos[j]+f3
                #     tp2 = test(middle2,middle_l,c)
                #     ts2 = test(middle2,middle_l,a)
                #     tt2 = test(middle2,middle_l,b)
                #     if(tp!=tp2):
                #         test(middle2,middle_l,pos[j],debug=True)
                #         test(middle2,middle_l,c,debug=True)

                #     assert (tp==tp2) and (tt==tt2) and (ts2==ts),f"\n{tp}\t{tp2}\t{ts}\t{ts2}\t{tt}\t{tt2}\n{pos[j]}\t{close_pos[j]}\t{far_pos[j]}\n{c}\t{a}\t{b}\n{middle2}\n{middle_l}\n{f1}\t{f2}\t{f3}"
                self.set_max_move(f,max_move,points[not_on_edge],angle,np.clip(close_dist/2.1-self.gamma/2,0,None),False)
                self.set_max_move(f,max_move,close,(angle+np.pi)%(np.pi*2),np.clip(close_dist/2.1-self.gamma/2,0,None),False)
                self.set_max_move(f,max_move,far,(angle+np.pi)%(np.pi*2),np.clip(np.linalg.norm(projs2-far_pos,axis=1)/1.05-self.gamma/2,0,None),False)
                # if i==21 and len(self.pos)>21:
                    # print("debug HSNC NvsIII b",self.pos[21],f[21],max_move[21],global_max,far,close)

                    
        return max_move

    def move(self,f,max_move):
        mag = np.linalg.norm(f,axis=1)
        # factor = np.where(mag>max_move,max_move/mag,np.ones(len(f)))
        to_change = mag>max_move
        # if(self.g.numberOfNodes()==380):
        #     print(f[249],"    ",f[0],max_move[249])
        f[to_change]*=np.expand_dims(max_move[to_change]/mag[to_change],-1)
        # if(self.g.numberOfNodes()==380):
        #     print(f[249],"    ",f[0])

        # f = f*np.expand_dims(factor,-1)
        self.pos+=f
        self.kdt = sspatial.KDTree(self.pos)
    def compute_faces(self):
        pass


if __name__ == "__main__":
    g = tlp.loadGraph("/tmp/testForce.tlpb.gz")
    vl = g.getLayoutProperty("viewLayout")
    vl.center()
    with open(f"/tmp/testforce.txt","wt") as f:
        f.write(f"{g.numberOfNodes()}\n")
        for i,n in enumerate(g.nodes()):
            p=vl[n]
            f.write(f"{i};{p.getX()};{p.getY()}\n")
        f.write(f"{g.numberOfEdges()}\n")
        for i,e in enumerate(g.edges()):
            src = g.nodePos(g.source(e))
            tgt = g.nodePos(g.target(e))
            f.write(f"{src};{tgt}\n")
    # r = tlp.computeBoundingRadius(g)
    # d = r[0].dist(r[1])
    # vl.scale(tlp.Vec3f(1/d,1/d,1.0))
    i = Impred(g,delta=2,gamma=5,nb_iter=1,zones=None,surrounding_edges={n:g.edges() for n in g.getNodes()},fixed_nodes=set())
    # for j in range(360):
    #     a = j*np.pi/180
    #     p2 = np.random.uniform(0,100)*np.stack((np.cos(a),np.sin(a)),-1)
    #     a2 = i.test_angle(p2[1],p2[0])
    #     b = np.abs(a-a2)<1e-2
    #     if(not b):
    #         i.test_angle(p2[1],p2[0],debug=True)
    #     assert b,(a,a2,j,a2*180/np.pi)
    # def test(a,b,c):
    #     return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])>0
    # for j in range(1,180):
    #     angle1 = np.random.uniform(0,np.pi*2)
    #     angle = angle1 + j*np.pi/180
    #     a = np.random.uniform(size=(2))
    #     b = a + np.stack((np.cos(angle1),np.sin(angle1)),-1)
    #     c = a + np.stack((np.cos(angle),np.sin(angle)),-1)*np.random.uniform()
    #     assert test(a,b,c),(a,b,c)
    # for j in range(181,360):
    #     angle1 = np.random.uniform(0,np.pi*2)
    #     angle = angle1 + j*np.pi/180
    #     a = np.random.uniform(size=(2))
    #     b = a + np.stack((np.cos(angle1),np.sin(angle1)),-1)
    #     c = a + np.stack((np.cos(angle),np.sin(angle)),-1)*np.random.uniform()
    #     assert not test(a,b,c),(a,b,c)
    pos = i.layout()
    for i,n in enumerate(g.nodes()):
        vl[n]=tlp.Vec3f(pos[i,0],pos[i,1],0.0)
    with open(f"/tmp/testforceo.txt","wt") as f:
        f.write(f"{g.numberOfNodes()}\n")
        for i,n in enumerate(g.nodes()):
            p=vl[n]
            f.write(f"{i};{p.getX()};{p.getY()}\n")
        f.write(f"{g.numberOfEdges()}\n")
        for i,e in enumerate(g.edges()):
            src = g.nodePos(g.source(e))
            tgt = g.nodePos(g.target(e))
            f.write(f"{src};{tgt}\n")

    # r = tlp.computeBoundingRadius(g)
    # d2 = r[0].dist(r[1])
    # vl.scale(tlp.Vec3f(d/d2,d/d2,1.0))
    
    tlp.saveGraph(g,"/tmp/testimpred.tlpb.gz")