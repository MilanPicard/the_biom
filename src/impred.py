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
        tlp.saveGraph(self.g,f"/tmp/start.tlpb.gz")

        self.delta = delta
        self.gamma = gamma
        self.nb_iter = nb_iter
        self.vl = g.getLayoutProperty("viewLayout")
        if(g.numberOfNodes()>100):
            print(self.vl[g.nodes()[98]],self.vl[g.nodes()[100]])

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
        if(len(pos)>100):
            print(self.pos[98],self.pos[100])
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
        assert np.all(np.isfinite(self.pos))

    def layout(self):
        self.compute_faces()
        global_maxes = np.linspace(3*self.delta,0,num=self.nb_iter)
        vl=self.g.getLayoutProperty("viewLayout")
        tlp.saveGraph(self.g,f"/tmp/iter_-1.tlpb.gz")
        assert np.all(np.isfinite(self.pos))
        for j,global_max in enumerate(global_maxes):
            temp = (self.nb_iter-j)/self.nb_iter
            temp=1
            self.iter(global_max,temp)
            for i,n in enumerate(self.g.nodes()):
                vl[n]=tlp.Vec3f(self.pos[i,0],self.pos[i,1],0.0)
            tlp.saveGraph(self.g,f"/tmp/iter_{j}.tlpb.gz")
            # if  j==0:
                # break
        
        vl.scale(tlp.Vec3f(self.scale,self.scale,1.0))
        vl.center(self.center)
        # return self.pos*self.scale + np.array([[self.center.getX(),self.center.getY()]])
        return self.pos + np.array([[self.center.getX(),self.center.getY()]])
    def iter(self,global_max,temp):
        t1 = time.perf_counter()
        f = self.compute_forces(temp)
        t2 = time.perf_counter()
        max_move=self.compute_max_move(global_max,f)
        t3 = time.perf_counter()
        # print(np.count_nonzero(max_move==global_max))
        # if len(self.pos)>21:
            # print("debug HSNC NvsIII",self.pos[21],f[21],max_move[21],global_max)
        self.move(f,max_move)
        t4 = time.perf_counter()
        print("iter",t2-t1,t3-t2,t4-t3)
    def compute_forces(self,temp):
        fg = self.gravity_force(temp)
        fn = self.node_node_rep(temp)
        fa = self.edge_attr()
        fe = self.edge_rep()
        r = fn+fa+fe+fg
        assert np.all(np.isfinite(fg))
        assert np.all(np.isfinite(fn))
        assert np.all(np.isfinite(fa))
        assert np.all(np.isfinite(fe))
        # print("mean rep",np.mean(np.linalg.norm(fn[self.is_movable],axis=1)),np.mean(np.linalg.norm(fe[self.is_movable],axis=1)),np.mean(np.linalg.norm(fa[self.is_movable],axis=1)),self.vl.averageEdgeLength())
        # r = fe
        # if(self.g.numberOfNodes()==380):
            # print(self.pos[249],fg[249],fn[249],fa[249],fe[249],r[249],"    ",self.pos[0],fg[0],fn[0],fa[0],fe[0],r[0])
        # for i in range(len(fn)):
        # print(np.dot(fn[0],r[0])/(np.linalg.norm(fn[0])*np.linalg.norm(r[0])),np.arctan2(fn[0,1],fn[0,0])*180/np.pi,np.arctan2(r[0,1],r[0,0])*180/np.pi)
        return r
    def gravity_force(self,temp):
        # c = self.pos[self.movable_index].mean(axis=0)
        # diff = np.expand_dims(c,0)-self.pos[self.movable_index]
        diff = -self.pos[self.movable_index]
        dist = np.linalg.norm(diff,axis=1)
        f = np.zeros_like(self.pos)
        f[self.movable_index] = diff*(np.expand_dims(dist/10,-1)**(1+(temp)))
        return f

    def node_node_rep(self,temp):
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
            # print("minDist",np.min(dist),np.mean(dist))
            mag = ((self.delta)/dist)**2/dist
            # print("aaaaaaa",self.g.numberOfEdges(),self.g.numberOfNodes(),len(pairs),mag.mean(),self.vl.averageEdgeLength())

            repulsions = np.expand_dims(mag,-1) *diff
            # repulsions = np.expand_dims((1/dist)**3,-1) *diff
            assert np.all(dist>0),(np.argwhere(dist<=0),pairs[np.argwhere(dist<=0)[0]],self.pos[np.unique(pairs[np.argwhere(dist<=0)[0]].flatten())])
            assert np.all(np.isfinite(repulsions))
            forces[:,0]+=np.bincount(pairs[first_movable,0],repulsions[first_movable,0],minlength=self.pos.shape[0])
            assert np.all(np.isfinite(forces))
            forces[:,1]+=np.bincount(pairs[first_movable,0],repulsions[first_movable,1],minlength=self.pos.shape[0])
            assert np.all(np.isfinite(forces))
            forces[:,0]+=np.bincount(pairs[second_movable,1],-repulsions[second_movable,0],minlength=self.pos.shape[0])
            assert np.all(np.isfinite(forces))
            forces[:,1]+=np.bincount(pairs[second_movable,1],-repulsions[second_movable,1],minlength=self.pos.shape[0])

        # if len(forces)==7:
            # print(forces,self.pos,pairs,repulsions)
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
        pos_src = []
        pos_tgt = []
        pos_point = []
        point_index = []
        src_index = []
        tgt_index = []
        for i,n in enumerate(self.movable_nodes):
            srcs = []
            tgts = []
            if n in self.surrounding_edges:
                for e in self.surrounding_edges[n]:
                    if not (n == self.g.source(e) or n == self.g.target(e)):
                        srcs.append(self.g.nodePos(self.g.source(e)))
                        tgts.append(self.g.nodePos(self.g.target(e)))
                pos_src.append(self.pos[srcs])
                pos_tgt.append(self.pos[tgts])
                j = self.g.nodePos(n)
                pos_point.append(np.repeat(np.expand_dims(self.pos[j],0),len(srcs),0))
                point_index.append(np.full(len(srcs),j))
                src_index.append(srcs)
                tgt_index.append(tgts)
        forces = np.zeros_like(self.pos)
        dist = None
        if len(pos_src)>0:
            pos_src = np.concatenate(pos_src)
            pos_tgt = np.concatenate(pos_tgt)
            pos_point = np.concatenate(pos_point)
            point_index = np.concatenate(point_index)
            src_index = np.concatenate(src_index)
            tgt_index = np.concatenate(tgt_index)
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
        # if np.any(point_index==176):
            # print("176     ",forces[176],pos_point[point_index==176], proj[point_index==176],proj_dist[point_index==176],pos_src[point_index==176],pos_tgt[point_index==176],close_enough[point_index==176])#,on_edge[close_enough[point_index==176]])
        # elif self.g.numberOfNodes()>10:
            # print(self.g.getStringProperty("viewLabel")[self.g.nodes()[176]],self.g.nodes()[176] in self.surrounding_edges,self.g.nodes()[176] in self.movable_nodes,np.unique(point_index ),len(np.unique(point_index )))
        # print(len(src_index),self.g.numberOfNodes(),dist)
        return forces

        
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
    def set_max_move(self,f,max_move,i,angle,distance,debug=False,times=None):
        t1 = time.perf_counter()
        cond1 = np.linalg.norm(f[i],axis=1)>distance
        i2 = i[cond1]
        do_something = np.any(cond1)
        t2 = time.perf_counter()
        times[0]+=t2-t1
        if do_something:
            # fa = np.arctan2(f[i,1],f[i,0])%(np.pi*2)
            # fa = np.array([self.test_angle(f[j,1],f[j,0]) for j in i])%(np.pi*2)
            fa = np.arctan2(f[i2,1],f[i2,0])
            fa[fa>np.pi]-=np.pi*2
            fa = fa%(np.pi*2)
            # assert np.max(np.abs(fa-fa2))<1e-3,np.max(np.abs(fa-fa2))
            t3 = time.perf_counter()
            times[1]+=t3-t2
            angle_diff = (fa-angle[cond1])%(np.pi*2)
            angle_diff = np.abs(np.where(angle_diff>np.pi,angle_diff-np.pi*2,angle_diff))
            cond = np.logical_or(angle_diff<np.pi/2,angle_diff>np.pi*3/2)
            do_something = np.any(cond)
            t4 = time.perf_counter()
            times[2]+=t4-t3

            if do_something:
                # assert np.all(angle_diff[angle_diff<np.pi/2]>=0)
                # assert np.all(angle_diff[angle_diff>np.pi*3/2]<np.pi*2)
                # to_change =i[cond]
                d = distance[cond1][cond]
                # assert(np.all(d/np.cos(angle_diff[cond])>=0)),(f[to_change],max_move[to_change],angle[cond]*180/np.pi,distance[cond],d/np.cos(angle_diff[cond]))
                # if(debug):
                #     # print(angle,fa,angle_diff,angle*180/np.pi,fa*180/np.pi,angle_diff*180/np.pi,cond)
                #     for k,j in enumerate(i):
                #         if j==249:
                #             fa2 = self.gateway.jvm.ocotillo.geometry.Geom2D.angle(self.gateway.jvm.ocotillo.geometry.Coordinates(f[j,0],f[j,1]))
                #             # assert np.abs(fa[k]-fa2)<1e-3,(fa[k],fa2,f[j],np.arctan2(f[j,1],f[j,0]))
                #             fca = self.gateway.jvm.Math.abs(self.gateway.jvm.ocotillo.geometry.Geom2D.normalizeRadiansAngle(fa2-angle[k]))
                #             assert np.abs(angle_diff[k]-fca)<1e-3,f"{angle_diff[k]},{fca}\n{(fa[k],angle[k])}\n{(fa2,angle[k])}\n{self.gateway.jvm.ocotillo.geometry.Geom2D.posNormalizeRadiansAngle(fa2-angle[k])}\n{fa[k]-angle[k]}"
                #             # print("azertyuiop",j,distance[k],cond[k],distance[k]/np.cos(angle_diff[k]),max_move[j])
                # max_move[to_change]=np.minimum(max_move[to_change],d/np.cos(angle_diff[cond]))
                t5 = time.perf_counter()
                times[3]+=t5-t4
                # m = pd.DataFrame.from_dict({"id":i2[cond],"value":d/np.cos(angle_diff[cond])}).groupby("id").min()
                i3 = i2[cond]
                order = np.argsort(i3)
                i3 = i3[order]
                v = (d/np.cos(angle_diff[cond]))[order]
                uniques = np.unique(i3,return_index=True)
                v=np.array([np.min(j) for j in np.split(v,uniques[1][1:])])
                # for j in range(len(uniques[0])):
                    # assert m.loc[uniques[0][j]]["value"]==v[j]
                max_move[uniques[0]]=np.minimum(max_move[uniques[0]],v)
                # index = m.index.to_numpy()
                # max_move[index]=np.minimum(max_move[index],m["value"].to_numpy())
                t6 = time.perf_counter()
                times[4]+=t6-t5
                # times[9]+=t3-t2

                # if(debug):
                    # print(angle,fa,angle_diff,angle*180/np.pi,fa*180/np.pi,angle_diff*180/np.pi,cond)
                    # for k,j in enumerate(i):
                    #     if j==249:
                    #         print("azertyuiop2",j,distance[k],cond[k],distance[k]/np.cos(angle_diff[k]),max_move[j])

        # times[8]+=t2-t1


    def compute_max_move(self,global_max,f):
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
        times = []
        times2 = np.zeros(10)
        tz2 = time.perf_counter()
        srcs= []
        tgts= []
        points = []
        for i,n in enumerate((self.g.nodes())):
            t1 = time.perf_counter()
            if not self.is_movable[i]:
                times.append(time.perf_counter()-t1)
                continue
            pos_srcs = []
            pos_tgts = []
            
            if n in self.surrounding_edges:
                tx1 = time.perf_counter()
                for e in self.surrounding_edges[n]:
                    
                    src = self.g.nodePos(self.g.source(e))
                    tgt = self.g.nodePos(self.g.target(e))
                    if i != src and i!=tgt:
                        srcs.append(src)
                        tgts.append(tgt)
                        points.append(i)
                tx2 = time.perf_counter()
                # times2[0]+=tx2-tx1
                if len(srcs)==0:
                    times.append(time.perf_counter()-t1)
                    continue
        if len(srcs)!=0:
            srcs = np.array(srcs)
            tgts = np.array(tgts)
            points = np.array(points)
            pos_srcs = self.pos[srcs]
            pos_tgts = self.pos[tgts]
            # pos = np.expand_dims(self.pos[i],0).repeat(len(pos_srcs),0)
            pos = self.pos[points]
            projs = self.get_projections(pos_srcs,pos_tgts,pos)
            dist_p = np.linalg.norm(projs-pos,axis=1)
            not_too_far = dist_p<=global_max*2
            tx3 = time.perf_counter()
            # times2[1]+=tx3-tx2
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

            tx4 = time.perf_counter()
            # times2[2]+=tx4-tx3
            if(np.any(on_edge)):
                # point_coord =self.gateway.jvm.ocotillo.geometry.Coordinates(pos[0,0],pos[0,1],0.0)
                ty1 = time.perf_counter()
                on_edge_projs = projs[on_edge]
                on_edge_pos = pos[on_edge]
                angle = np.arctan2(on_edge_projs[:,1]-on_edge_pos[:,1],on_edge_projs[:,0]-on_edge_pos[:,0])%(2*np.pi)
                on_edge_dist = dist_p[on_edge]
                ty2 = time.perf_counter()
                # times2[5]+=ty2-ty1

                # for j in range(len(on_edge_pos)):
                # #     src_coord =self.gateway.jvm.ocotillo.geometry.Coordinates(pos_srcs[on_edge][j,0],pos_srcs[on_edge][j,1],0.0)
                # #     tgt_coord =self.gateway.jvm.ocotillo.geometry.Coordinates(pos_tgts[on_edge][j,0],pos_tgts[on_edge][j,1],0.0)
                # #     proj_coord =self.gateway.jvm.ocotillo.geometry.Geom2D.pointOnLineProjection(point_coord,src_coord,tgt_coord)
                # #     # assert np.linalg.norm(on_edge_projs[j]-np.array([proj_coord.x(),proj_coord.y()]))<1e-6#,(on_edge_projs[j],proj_coord.x(),proj_coord.y(),self.get_projections(pos_srcs[on_edge][j:j+1],pos_tgts[on_edge][j:j+1],pos[on_edge][j:j+1], debug=True))
                # #     pn = proj_coord.minus(point_coord)
                # #     angle2 = self.gateway.jvm.ocotillo.geometry.Geom2D.angle(pn)
                # #     assert np.abs(angle[j]-angle2)<1e-3,(angle,angle2)
                #     middle = (on_edge_pos[j]+on_edge_projs[j])/2
                #     middle_l = middle + np.stack((np.cos(angle[j]+np.pi/2),np.sin(angle[j]+np.pi/2)))
                #     # t1 = test(pos_srcs[on_edge][j],pos_tgts[on_edge][j],[j])
                #     tp = test(middle,middle_l,on_edge_pos[j])
                #     ts = test(middle,middle_l,pos_srcs[on_edge][j])
                #     tt = test(middle,middle_l,pos_tgts[on_edge][j])
                #     self.set_max_move(f,max_move,np.full((1,),i),angle[j:j+1],on_edge_dist[j:j+1]/2.1)
                #     self.set_max_move(f,max_move,srcs[on_edge][j:j+1],angle[j:j+1]+np.pi,on_edge_dist[j:j+1]/2.1)
                #     self.set_max_move(f,max_move,tgts[on_edge][j:j+1],angle[j:j+1]+np.pi,on_edge_dist[j:j+1]/2.1)
                #     mag = np.linalg.norm(f[i],keepdims=True)
                #     if mag[0]>max_move[i]:
                #         f1 = f[i]*max_move[i]/mag
                #     else:
                #         f1 = f[i]
                #     c = self.pos[i]+f1
                #     mag = np.linalg.norm(f[srcs[on_edge][j]],keepdims=True)
                #     if mag[0]>max_move[srcs[on_edge][j]]:
                #         f2 = f[srcs[on_edge][j]]*max_move[srcs[on_edge][j]]/mag
                #     else:
                #         f2 = f[srcs[on_edge][j]]
                #     a = self.pos[srcs[on_edge][j]]+f2
                    
                #     mag = np.linalg.norm(f[tgts[on_edge][j]],keepdims=True)
                #     if mag[0]>max_move[tgts[on_edge][j]]:
                #         f3 = f[tgts[on_edge][j]]*max_move[tgts[on_edge][j]]/mag
                #     else:
                #         f3 = f[tgts[on_edge][j]]
                #     b = self.pos[tgts[on_edge][j]]+f3
                #     tp2 = test(middle,middle_l,c)
                #     ts2 = test(middle,middle_l,a)
                #     tt2 = test(middle,middle_l,b)
                #     if(tp!=tp2):
                #         test(middle,middle_l,on_edge_pos[j],debug=True)
                #         test(middle,middle_l,c,debug=True)

                #     assert (tp==tp2) and (tt==tt2) and (ts2==ts),f"\n{tp}\t{tp2}\t{ts}\t{ts2}\t{tt}\t{tt2}\n{on_edge_pos[j]}\t{pos_srcs[on_edge][j]}\t{pos_tgts[on_edge][j]}\n{c}\t{a}\t{b}\n{on_edge_projs[j]}\n{middle}\n{middle_l}"
                    
                    # t2 = test(a,b,c)
                    # if(t1!=t2):
                    #     self.set_max_move(f,max_move,np.full((1,),i),angle[j:j+1],on_edge_dist[j:j+1],debug=True)
                    #     self.set_max_move(f,max_move,srcs[on_edge][j:j+1],angle[j:j+1]+np.pi,on_edge_dist[j:j+1]/2,True)
                    #     self.set_max_move(f,max_move,tgts[on_edge][j:j+1],angle[j:j+1]+np.pi,on_edge_dist[j:j+1]/2,True)

                    # assert (t1==t2)==(np.dot(a-b,self.pos[srcs[on_edge][j]]-self.pos[tgts[on_edge][j]])>=0),(t1,t2,self.pos[i],self.pos[srcs[on_edge][j]],self.pos[tgts[on_edge][j]],
                                #    f1,f2,f3,
                                #    a,b,c,
                                #    on_edge_projs[j])

                self.set_max_move(f,max_move,points[on_edge],angle,np.clip(on_edge_dist/2.1-self.gamma/2,0,None),debug=False,times=times2)
                self.set_max_move(f,max_move,srcs[on_edge],angle+np.pi,np.clip(on_edge_dist/2.1-self.gamma/2,0,None),debug=False,times=times2)
                self.set_max_move(f,max_move,tgts[on_edge],angle+np.pi,np.clip(on_edge_dist/2.1-self.gamma/2,0,None),debug=False,times=times2)
                # if i==21 and len(self.pos)>21:
                    # print("debug HSNC NvsIII a",self.pos[21],f[21],max_move[21],global_max)

                ty3 = time.perf_counter()
                # times2[6]+=ty3-ty2
                tx5 = time.perf_counter()
                # times2[3]+=tx5-tx4
            if(not np.all(on_edge)):
                ty4 = time.perf_counter()
                
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
                ty5 = time.perf_counter()
                # times2[7]+=ty5-ty4
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
                self.set_max_move(f,max_move,points[not_on_edge],angle,np.clip(close_dist/2.1-self.gamma/2,0,None),False,times=times2)
                self.set_max_move(f,max_move,close,(angle+np.pi)%(np.pi*2),np.clip(close_dist/2.1-self.gamma/2,0,None),False,times=times2)
                self.set_max_move(f,max_move,far,(angle+np.pi)%(np.pi*2),np.clip(np.linalg.norm(projs2-far_pos,axis=1)/1.05-self.gamma/2,0,None),False,times=times2)
                # if i==21 and len(self.pos)>21:
                    # print("debug HSNC NvsIII b",self.pos[21],f[21],max_move[21],global_max,far,close)

                ty6 = time.perf_counter()
                # times2[8]+=ty6-ty5
                tx6 = time.perf_counter()
                # times2[4]+=tx6-tx5
                    
            # times.append(time.perf_counter()-t1)
        # print(pd.DataFrame.from_dict({"i":np.arange(len(times)),"t":np.array(times)}).sort_values("t",ascending=False).head(),times2)
        tz3 = time.perf_counter()
        # print("times2",times2,tz2-tz1,tz3-tz2)
        return max_move

    def move(self,f,max_move):
        assert np.all(np.isfinite(self.pos))
        assert np.all(np.isfinite(f))
        assert np.all(np.isfinite(max_move))
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