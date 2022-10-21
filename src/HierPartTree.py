import numpy as np
import pandas as pd
import time as t

def heterogeneity(vertex):
    return HPTree.dist[vertex].loc[vertex].sum().sum() / (len(vertex) * 2)

def combine_best_h(t1,t2):
    new_best_h = {}
    for i in t1.best_h.keys():
        for j in t2.best_h.keys():
            h = t1.best_h[i][0] + t2.best_h[j][0]
            if i+j in new_best_h:
                if h < new_best_h[(i+j)][0]:
                    new_best_h[(i + j)] = (h,t1.best_h[i][1] + t2.best_h[j][1])
            else:
                new_best_h[(i + j)] = (h,t1.best_h[i][1] + t2.best_h[j][1])
    return new_best_h

def clearHPT():
    HPTree.dict_regions = {}
    HPTree.dist = None
    HPTree.identical_nodes = 0


class HPTree():
    dist = None
    dict_regions = {}
    identical_nodes = 0
    def __init__(self, vertex, neighbors, edges, h, center=None):
        self.vertex = vertex
        self.neighbors = neighbors
        self.h = h
        self.edges = {}
        self.center = center
        self.best_h = {1:(h,[self])}
        self.children = []
        if type(edges) is list:
            self.create_arch(edges)
        else:
            self.edges = edges
        if self.center is None:
            self.center = {}
            for e in self.edges.keys():
                self.center[e] = (len(self.edges[e][e[0]]),len(self.edges[e][e[1]]))
        if len(self.vertex) > 1:
            self.middle_edge = max(self.center.keys(),key=(lambda e: min(self.center[e])))

    def copy(self):
        return HPTree(self.vertex.copy(),self.neighbors.copy(),self.edges.copy(),self.h,self.center.copy())

    def h_after_cut(self,e):
        #t1 = set(self.vertex).intersection(self.edges[e][e[0]])
        #t2 = set(self.vertex).intersection(self.edges[e][e[1]])
        h1, h2 = heterogeneity(self.edges[e][e[0]]), heterogeneity(self.edges[e][e[1]])#h1, h2 = heterogeneity(list(t1)), heterogeneity(list(t2))
        return h1, h2

    def create_arch(self,T):
        for e in T:
            e = (min(e),max(e))
            self.edges[e] = {e[0]:[e[0]],e[1]:[e[1]]}
            self.collect(e,e[0],e[1])
            self.collect(e,e[1],e[0])
            self.spread(self.edges[e][e[0]],e[1],e[0])
            self.spread(self.edges[e][e[1]],e[0],e[1])

    def collect(self,e,To,From):
        for v in self.neighbors[To]:
            if v != From:
                new_e = (min(To,v),max(To,v))
                if new_e in self.edges:
                    self.edges[e][To] += self.edges[new_e][v]

    def spread(self,val,To,From):
        for v in self.neighbors[To]:
            if v != From:
                new_e = (min(To,v),max(To,v))
                if new_e in self.edges:
                    self.edges[new_e][To] += val
                    self.spread(val,v,To)

    def compute_best_h(self,w,k):
        dicts_h = []
        for i in range(w):
            dicts_h.append(combine_best_h(self.children[i * 2], self.children[i * 2 + 1]))
        l = 2
        ok = True
        while ok and l <= k:
            best_d = None
            for d in dicts_h:
                if l in d and (best_d is None or d[l][0] < best_d[l][0]):
                    best_d = d
            if best_d is None:
                ok = False
            else:
                self.best_h[l] = best_d[l]
                l += 1

    def find_regions(self, k, limit_depth, w_list, sc_list, n_best, index_w=0,index_sc=0, depth=0,talk=False):
        if str(self.vertex) in HPTree.dict_regions:
            HPTree.identical_nodes += 1
        else:
            if talk:
                if depth >= 2:
                    talk = False
                print()
                print('Node of size '+str(len(self.vertex))+ ' and depth = '+str(depth))
            HPTree.dict_regions[str(self.vertex)] = self
            if len(self.vertex) > 1 and depth < limit_depth:
                if index_w < len(w_list)-1 and depth == w_list[index_w+1][0]:
                    index_w += 1
                w = w_list[index_w][1]
                if index_sc < len(sc_list)-1 and depth == sc_list[index_sc+1][0]:
                    index_sc += 1
                sc = sc_list[index_sc][1]
                if depth == limit_depth-1:
                    w=1
                new_trees = self.tree_partitioning(w,sc,n_best,talk=talk)
                for tree in new_trees:
                    tree.find_regions(k, limit_depth, w_list, sc_list, n_best, index_w, index_sc, depth=depth+1,talk=talk)
                    self.children.append(HPTree.dict_regions[str(tree.vertex)])
                self.compute_best_h(w,k)




    def tree_partitioning(self,w,sc,n_best,talk=False):
        if talk:
            t1 = t.time()
        e,h,h1,h2 = self.find_best_cut(w,sc,n_best=n_best)
        if talk:
            t2 = t.time() - t1
            print('time to find ' + str(w) + ' best cuts with sc = '+str(sc) + ' : ' + str(t2))
        if talk:
            t1 = t.time()
        new_trees = []
        for i in range(w):
            neighbors = self.neighbors.copy()
            #vertex1 = list(set(self.vertex).intersection(self.edges[e[i]][e[i][0]]))
            #vertex2 = list(set(self.vertex).intersection(self.edges[e[i]][e[i][1]]))
            vertex1 = self.edges[e[i]][e[i][0]]
            vertex2 = self.edges[e[i]][e[i][1]]
            center1, center2 = {}, {}
            edges1, edges2 = {}, {}
            for edge in self.edges.keys():#for edge in self.center.keys():
                if edge != e[i]:
                    if edge[0] in vertex1:
                        if e[i][0] in self.edges[edge][edge[0]]:
                            edges1[edge] = {edge[0]:[v for v in self.edges[edge][edge[0]] if v in vertex1],edge[1]:self.edges[edge][edge[1]].copy()}
                            #center1[edge] = (self.center[edge][0]-len(vertex2),self.center[edge][1])
                        else:
                            edges1[edge] = {edge[0]: self.edges[edge][edge[0]].copy(),
                                            edge[1]: [v for v in self.edges[edge][edge[1]] if v in vertex1]}
                            #center1[edge] = (self.center[edge][0],self.center[edge][1]-len(vertex2))
                        center1[edge] = (len(edges1[edge][edge[0]]),len(edges1[edge][edge[1]]))
                    else:
                        if e[i][0] in self.edges[edge][edge[0]]:
                            edges2[edge] = {edge[0]: [v for v in self.edges[edge][edge[0]] if v in vertex2],
                                            edge[1]: self.edges[edge][edge[1]].copy()}
                            #center2[edge] = (self.center[edge][0]-len(vertex1),self.center[edge][1])
                        else:
                            edges2[edge] = {edge[0]: self.edges[edge][edge[0]].copy(),
                                            edge[1]: [v for v in self.edges[edge][edge[1]] if v in vertex2]}
                            #center2[edge] = (self.center[edge][0],self.center[edge][1]-len(vertex1))
                        center2[edge] = (len(edges2[edge][edge[0]]), len(edges2[edge][edge[1]]))
            neighbors[e[i][0]] = [v for v in neighbors[e[i][0]] if v!=e[i][1]]
            neighbors[e[i][1]] = [v for v in neighbors[e[i][1]] if v!=e[i][0]]
            if str(vertex1) in HPTree.dict_regions:
                new_trees += [HPTree.dict_regions[str(vertex1)]]
            else:
                new_trees += [HPTree(vertex1,neighbors,edges1,h1[i], center=center1)]
            if str(vertex2) in HPTree.dict_regions:
                new_trees += [HPTree.dict_regions[str(vertex2)]]
            else:
                new_trees += [HPTree(vertex2, neighbors, edges2, h2[i], center=center2)]
        if talk:
            t2 = t.time() - t1
            print('time to create ' + str(w) + ' new trees : ' + str(t2))
        return new_trees

    def find_best_cut(self, w, sc, n_best=5):
        if w==1:
            n_best=1
        best_h1, best_h2 = [0.0 for _ in range(n_best)], [0.0 for _ in range(n_best)]
        best_h1[0], best_h2[0] = self.h_after_cut(self.middle_edge)
        best_cut = [self.middle_edge for _ in range(n_best)]
        best_h = [float('inf') for _ in range(n_best)]
        best_h[0] = best_h1[0] + best_h2[0]
        max_best_h = max(best_h)
        index_max_h = best_h.index(max_best_h)
        calculated = [self.middle_edge]
        i_since_upgrade = 0
        L = {}
        to_expand = self.middle_edge
        expand = True
        while i_since_upgrade < sc and expand:
            best_neighbor_h = float('inf')
            new_best = False
            for b1 in to_expand:
                for b2 in self.neighbors[b1]:
                    e = (min(b1,b2),max(b1,b2))
                    if e not in calculated:
                        calculated += [e]
                        h1, h2 = self.h_after_cut(e)
                        h = h1 + h2
                        L[e] = max(h1,h2)
                        if h < max_best_h and h < best_neighbor_h:
                            best_neighbor_h = h
                            new_best = True
                            best_h[index_max_h], best_cut[index_max_h] = h, e
                            best_h1[index_max_h], best_h2[index_max_h] = h1, h2
                            # max_best_h = max(best_h)
                            # index_max_h = best_h.index(max_best_h)
                            # i_since_upgrade = 0
            i_since_upgrade += 1
            if new_best:
                max_best_h = max(best_h)
                index_max_h = best_h.index(max_best_h)
                i_since_upgrade = 0
            if len(L) > 0:
                to_expand = min(L,key=L.get)
                del L[to_expand]
            else:
                expand = False
        w_best_h1, w_best_h2 = [0.0 for _ in range(w)], [0.0 for _ in range(w)]
        w_best_cut = [0.0 for _ in range(w)]
        w_best_h = [0.0 for _ in range(w)]
        min_best_h = min(best_h)
        i_min = best_h.index(min_best_h)
        w_best_h1[0], w_best_h2[0], w_best_cut[0], w_best_h[0] = best_h1[i_min], best_h2[i_min], best_cut[i_min], best_h[i_min]
        for i in range(1,w):
            max_diff = -1
            i_max_diff = 0
            for k in range(n_best):
                min_diff = float('inf')
                for j in range(i):
                    diff = abs(min(self.center[best_cut[k]])-min(self.center[w_best_cut[j]]))
                    min_diff = min(min_diff,diff)
                if min_diff > max_diff:
                    max_diff = min_diff
                    i_max_diff = k
            w_best_h1[i], w_best_h2[i], w_best_cut[i], w_best_h[i] = best_h1[i_max_diff], best_h2[i_max_diff], best_cut[i_max_diff], \
                                                                     best_h[i_max_diff]

        return w_best_cut, w_best_h, w_best_h1, w_best_h2
