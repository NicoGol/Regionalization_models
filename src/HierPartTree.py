import numpy as np
import pandas as pd
import time as t

def heterogeneity(vertex):
    return HPTree.dist[vertex].loc[vertex].sum().sum() / (len(vertex) * 2)


class HPTree():
    dist = None
    n_partitions = 0
    n_cuts = 0
    he = [[]]
    regions = [[]]
    dict_regions = {}
    dict_he = {}
    index_node = 0
    arch = {}
    def __init__(self, vertex, neighbors, edges, h, middle_edge, center=None):
        self.vertex = vertex
        self.neighbors = neighbors
        self.h = h
        self.edges = {}
        self.center = center
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
        # self.middle_edge = middle_edge
        # if self.middle_edge is None:
        #     middle_edge, best = None, 0
        #     for edge in self.edges.keys():
        #         middle = min(len(self.edges[edge][min(edge)]),len(self.edges[edge][max(edge)]))
        #         if middle > best:
        #             middle_edge, best = edge, middle
        #     self.middle_edge = middle_edge

    def h_after_cut(self,e):
        t1 = list(set(self.vertex).intersection(self.edges[e][e[0]]))
        t2 = list(set(self.vertex).intersection(self.edges[e][e[1]]))
        #return heterogeneity(self.edges[e][e[0]]), heterogeneity(self.edges[e][e[1]])
        h1, h2 = heterogeneity(t1), heterogeneity(t2)
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

    def find_regions(self, limit_vertex, limit_depth, w, n_sc, subtree, n_best, actual_w=0,actual_sc=0,actual_mins=0, depth=0, split_tree=False, index=0, min_size=0.0,talk=False):
        # if split_tree and depth == 0:
        #     Tree.regions = [[] for _ in range(n_split[actual_split][1])]
        #     Tree.he = [[] for _ in range(n_split[actual_split][1])]
        #Tree.regions[subtree] += [self.vertex]
        if talk:
            if depth >= 2:
                talk = False
            print()
            print('Node '+str(index)+ ' of size '+str(len(self.vertex))+ ' and depth = '+str(depth))
        HPTree.dict_regions[index] = self.vertex
        #Tree.he[subtree] += [self.h]
        HPTree.dict_he[index] = self.h
        if len(self.vertex) > limit_vertex * len(HPTree.dist) and depth < limit_depth:
            if actual_w < len(w)-1 and len(self.vertex) <= w[actual_w+1][0]*len(HPTree.dist):
                actual_w += 1
            if actual_sc < len(n_sc)-1 and len(self.vertex) <= n_sc[actual_sc+1][0]*len(HPTree.dist):
                actual_sc += 1
            if actual_mins < len(min_size)-1 and len(self.vertex) <= min_size[actual_mins+1][0]*len(HPTree.dist):
                actual_mins += 1
            if depth == limit_depth-1:
                new_trees = self.tree_partitioning(1, n_sc[actual_sc][1], min_size[actual_mins][1], n_best,
                                                   talk=talk)
            else:
                new_trees = self.tree_partitioning(w[actual_w][1],n_sc[actual_sc][1],min_size[actual_mins][1],n_best,talk=talk)
            new_indices = list(range(HPTree.index_node+1,HPTree.index_node+1+len(new_trees)))
            if talk:
                print('children of node '+str(index)+' are : '+str(new_indices))
            HPTree.arch[index] = new_indices
            HPTree.index_node += len(new_trees)
            HPTree.n_partitions += 1
            for i,tree in enumerate(new_trees):
                if split_tree and depth == 0:
                    tree.find_regions(limit_vertex, limit_depth, w, n_sc, int(i/2), n_best, actual_w, actual_sc, actual_mins, depth=depth+1, index=new_indices[i], min_size=min_size,talk=talk)
                else:
                    tree.find_regions(limit_vertex, limit_depth, w, n_sc, subtree, n_best, actual_w, actual_sc, actual_mins, depth=depth+1, index=new_indices[i], min_size=min_size,talk=talk)

    def tree_partitioning(self,z,sc,min_size,n_best,talk=False):
        if talk:
            t1 = t.time()
        e,h,h1,h2 = self.find_best_cut(z,sc,min_size,n_best=n_best)
        if talk:
            t2 = t.time() - t1
            print('time to find ' + str(z) + ' best cuts with sc = '+str(sc) + ' : ' + str(t2))
        if talk:
            t1 = t.time()
        new_trees = []
        HPTree.n_cuts += z
        for i in range(z):
            #tree_edges = self.edges.copy()
            neighbors = self.neighbors.copy()
            #vertex1, vertex2 = tree_edges[e[i]][e[i][0]], tree_edges[e[i]][e[i][1]]
            vertex1 = list(set(self.vertex).intersection(self.edges[e[i]][e[i][0]]))
            vertex2 = list(set(self.vertex).intersection(self.edges[e[i]][e[i][1]]))
            center1, center2 = {}, {}
            for edge in self.center.keys():
                if edge != e[i]:
                    if edge[0] in vertex1:
                        if e[i][0] in self.edges[edge][edge[0]]:
                            center1[edge] = (self.center[edge][0]-len(vertex2),self.center[edge][1])
                        else:
                            center1[edge] = (self.center[edge][0],self.center[edge][1]-len(vertex2))
                    else:
                        if e[i][0] in self.edges[edge][edge[0]]:
                            center2[edge] = (self.center[edge][0]-len(vertex1),self.center[edge][1])
                        else:
                            center2[edge] = (self.center[edge][0],self.center[edge][1]-len(vertex1))
            #edges1, edges2 = {}, {}
            neighbors[e[i][0]] = [v for v in neighbors[e[i][0]] if v!=e[i][1]]
            neighbors[e[i][1]] = [v for v in neighbors[e[i][1]] if v!=e[i][0]]
            # del tree_edges[e[i]]
            # middle_edge1, middle_edge2 = None, None
            # best_middle1, best_middle2 = 0, 0
            # for edge in tree_edges.keys():
            #     if edge[0] in vertex1:
            #         edges1[edge] = {edge[0]:[v for v in tree_edges[edge][edge[0]] if v in vertex1],
            #                         edge[1]:[v for v in tree_edges[edge][edge[1]] if v in vertex1]}
            #         middle = min(len(edges1[edge][edge[0]]),len(edges1[edge][edge[1]]))
            #         if middle > best_middle1:
            #             best_middle1 = middle
            #             middle_edge1 = edge
            #     else:
            #         edges2[edge] = {edge[0]: [v for v in tree_edges[edge][edge[0]] if v in vertex2],
            #                         edge[1]: [v for v in tree_edges[edge][edge[1]] if v in vertex2]}
            #         middle = min(len(edges2[edge][edge[0]]), len(edges2[edge][edge[1]]))
            #         if middle > best_middle2:
            #             best_middle2 = middle
            #             middle_edge2 = edge
            # new_trees += [Tree(vertex1,neighbors,edges1,h1[i],middle_edge1)]
            # new_trees += [Tree(vertex2, neighbors, edges2, h2[i], middle_edge2)]
            new_trees += [HPTree(vertex1,neighbors,self.edges,h1[i],None, center=center1)]
            new_trees += [HPTree(vertex2, neighbors, self.edges, h2[i], None, center=center2)]
        if talk:
            t2 = t.time() - t1
            print('time to create ' + str(z) + ' new trees : ' + str(t2))
        return new_trees

    def find_best_cut(self, z, sc, min_size, n_best=10):
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
        while i_since_upgrade < sc and expand:#len(calculated) < len(self.edges):
            best_neighbor_h = float('inf')
            new_best = False
            for v1 in to_expand:
                for v2 in self.neighbors[v1]:
                    e = (min(v1,v2),max(v1,v2))
                    if e not in calculated:
                        calculated += [e]
                        h1, h2 = self.h_after_cut(e)
                        h = h1 + h2
                        if min(self.center[e]) > len(self.vertex) * min_size:
                            L[e] = max(h1,h2)
                        if h < max_best_h and h < best_neighbor_h:
                            best_neighbor_h = h
                            best_h[index_max_h], best_cut[index_max_h] = h, e
                            best_h1[index_max_h], best_h2[index_max_h] = h1, h2
                            new_best = True
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
        z_best_h1, z_best_h2 = [0.0 for _ in range(z)], [0.0 for _ in range(z)]
        z_best_cut = [0.0 for _ in range(z)]
        z_best_h = [0.0 for _ in range(z)]
        min_best_h = min(best_h)
        i_min = best_h.index(min_best_h)
        z_best_h1[0], z_best_h2[0], z_best_cut[0], z_best_h[0] = best_h1[i_min], best_h2[i_min], best_cut[i_min], best_h[i_min]
        for i in range(1,z):
            max_diff = -1
            i_max_diff = 0
            for k in range(n_best):
                min_diff = float('inf')
                for j in range(i):
                    diff = abs(min(self.center[best_cut[k]])-min(self.center[z_best_cut[j]]))
                    min_diff = min(min_diff,diff)
                if min_diff > max_diff:
                    max_diff = min_diff
                    i_max_diff = k
            z_best_h1[i], z_best_h2[i], z_best_cut[i], z_best_h[i] = best_h1[i_max_diff], best_h2[i_max_diff], best_cut[i_max_diff], \
                                                                     best_h[i_max_diff]

        return z_best_cut, z_best_h, z_best_h1, z_best_h2

    def regions_to_onehot(self):
        # R = [pd.DataFrame(0,index=self.vertex,columns=range(len(Tree.regions[i]))) for i in range(len(Tree.regions))]
        # for j,regions in enumerate(Tree.regions):
        #     for i in range(len(regions)):
        #         for el in regions[i]:
        #             R[j].loc[el,i] = 1
        R = pd.DataFrame(0, index=self.vertex, columns=HPTree.dict_regions.keys())
        for i in HPTree.dict_regions.keys():
            for el in HPTree.dict_regions[i]:
                R.loc[el,i] = 1
        return R

    def clear(self):
        HPTree.regions = [[]]
        HPTree.n_cuts = 0
        HPTree.n_partitions = 0
        HPTree.he = [[]]
        HPTree.index_node = 0
        HPTree.arch = {}
        HPTree.dict_regions = {}
        HPTree.dict_he = {}