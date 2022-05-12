import numpy as np
import pandas as pd
import time as t
from sklearn.preprocessing import OneHotEncoder
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp

class Graph:
    def __init__(self, df, cols, full_order=None, dist_deg=2, talk=False):
        self.vertex = df.index.to_list()
        self.talk = talk
        self.times = {}
        if talk:
            print('start first_order')
        t1 = t.time()
        self.FirstO_E, self.contiguity_m, self.clusters = None, None, None
        self.first_order_edges(df, cols, dist_deg)
        t2 = t.time()
        self.times['FiO edges'] = t2-t1
        if talk:
            print('end first_order ' + str(self.times['FiO edges']))
        self.FullO_E, self.dist = {}, None
        if full_order is None:
            if talk:
                print('start full_order')
            t1 = t.time()
            self.full_order_edges(df, cols, dist_deg)
            t2 = t.time()
            self.times['FuO edges'] = t2 - t1
            if talk:
                print('end full_order ' + str(self.times['FuO edges']))
                print('start sort full order')
            t1 = t.time()
            self.sort_FullO_E()
            t2 = t.time()
            self.times['Sort FuO edges'] = t2 - t1
            if talk:
                print('stop sort full order ' + str(self.times['Sort FuO edges']))
        else:
            self.FullO_E, self.dist = full_order
        Tree.dist = self.dist
        self.dist_C = self.dist.copy()
        self.MST = None

    def first_order_edges(self,df, cols, dist_deg=2):
        First_order_edges = {}
        contiguity_matrix = pd.DataFrame(0, index=df.index, columns=df.index)
        clusters = {}
        for index, row in df.iterrows():
            e = {}
            for n in row.neighbors:
                dist = 0
                for col in cols:
                    dist += (row[col] - df.loc[n, col]) ** dist_deg
                # dist = np.linalg.norm(np.array(row[cols]) - np.array(df.loc[n, cols]))
                e[n] = {index: dist}
                contiguity_matrix.loc[index, n] = 1
                if index < n:
                    First_order_edges[(index, n)] = dist
            clusters[index] = Cluster(index, e)
        self.FirstO_E, self.contiguity_m, self.clusters = First_order_edges, contiguity_matrix, clusters

    def full_order_edges(self,df, cols, dist_deg=2):
        Full_order_edges = {}
        indexes = df.index
        edge_matrix = pd.DataFrame(0, index=indexes, columns=indexes)
        for i in range(len(indexes) - 1):
            index1 = indexes[i]
            for j in range(i + 1, len(indexes)):
                index2 = indexes[j]
                dist = 0
                for col in cols:
                    dist += (df.loc[index1, col] - df.loc[index2, col]) ** dist_deg
                edge_matrix.loc[index1, index2] = dist
                edge_matrix.loc[index2, index1] = dist
                if index1 < index2:
                    Full_order_edges[(index1, index2)] = dist
                else:
                    Full_order_edges[(index2, index1)] = dist
        self.FullO_E, self.dist = Full_order_edges, edge_matrix

    def sort_FirstO_E(self):
        self.FirstO_E = {k:v for k,v in sorted(self.FirstO_E.items(), key=lambda item: item[1])}

    def sort_FullO_E(self):
        self.FullO_E = {k:v for k,v in sorted(self.FullO_E.items(), key=lambda item: item[1])}

    def merge(self,l,m):
        self.clusters[m].go_off()
        self.clusters[l].grow(self.clusters[m])
        for v in self.clusters[m].get_v():
            self.clusters[v] = self.clusters[l]
        self.contiguity_m[l] = self.contiguity_m[[l,m]].max(axis=1)
        self.contiguity_m.loc[l] = self.contiguity_m.loc[[l, m]].max(axis=0)
        self.dist_C[l] = self.dist_C[[l, m]].max(axis=1)
        self.dist_C.loc[l] = self.dist_C.loc[[l, m]].max(axis=0)


    def find_MST_full_order_CL(self):
        i = 0
        T = {}
        edges = {}
        if self.talk:
            print('start MST')
        t1 = t.time()
        for (u,v) in self.FullO_E.keys():
            l, m = self.clusters[u].get_name(), self.clusters[v].get_name()
            if l != m and self.contiguity_m[l][m] != 0 and self.FullO_E[(u,v)] >= self.dist_C[l][m]:
                e,cost = self.clusters[u].find_shortest_edge(self.clusters[v])
                T[e] = cost
                self.merge(l,m) # changing clusters, changing contiguity and changing dist_C
                if e[0] in edges:
                    edges[e[0]] += [e[1]]
                else:
                    edges[e[0]] = [e[1]]
                if e[1] in edges:
                    edges[e[1]] += [e[0]]
                else:
                    edges[e[1]] = [e[0]]
                i += 1
        t2 = t.time()
        self.times['MST'] = t2-t1
        if self.talk:
            print('stop MST ' + str(self.times['MST']))
            print('start Tree arch')
        t1 = t.time()
        self.MST = Tree(self.vertex,edges,list(T.keys()),0,None)
        self.MST.h = heterogeneity(self.vertex)
        t2 = t.time()
        self.times['Tree arch'] = t2 - t1
        if self.talk:
            print('stop Tree arch ' + str(self.times['Tree arch']))
        return T

    def select_k_regions(self,k):
        if self.talk:
            print('start selecting regions')
        t1 = t.time()
        R = self.MST.regions_to_onehot()
        costs = np.array(self.MST.he)
        with gp.Env(empty=True) as env:
            env.setParam('OutputFlag', 0)
            env.start()
            with gp.Model(env=env) as model:
                n_regions = len(self.MST.regions)
                n_points = len(self.MST.vertex)
                x = model.addMVar(shape=n_regions, vtype=GRB.BINARY)
                model.setObjective(costs @ x)
                model.addConstr(x.sum() == k)
                sparse_R = sp.csr_matrix(R.values)
                val = np.array([1 for _ in range(n_points)])
                model.addConstr(sparse_R @ x == val)
                model.optimize()
                t2 = t.time()
                self.times['select regions'] = t2-t1
                if self.talk:
                    print('end selecting regions ' + str(self.times['select regions']))
                regions = [self.MST.regions[i] for i,val in enumerate(x.X) if val > 0.0]
                return model.objVal, regions

    def best_k_regions(self,k,limit_vertex, z, sc):
        T = self.find_MST_full_order_CL()
        if self.talk:
            print('start finding regions')
        t1 = t.time()
        self.MST.find_regions(limit_vertex, z, sc)
        t2 = t.time()
        self.times['column generation'] = t2-t1
        if self.talk:
            print('end finding regions ' + str(self.times['column generation']))
        res = self.select_k_regions(k)
        return res








class Cluster:
    def __init__(self,index,edges):
        self.name = index
        self.v = [index]
        self.e = edges
        self.on = True

    def get_name(self):
        return self.name

    def get_e(self):
        return self.e

    def get_v(self):
        return self.v

    def find_shortest_edge(self,cluster):
        best_u = self.v[0]
        best_dist = 1
        best_v = cluster.get_v()[0]
        for v in cluster.get_v():
            if v in self.e:
                for u in self.e[v].keys():
                    if self.e[v][u] < best_dist:
                        best_u = u
                        best_v = v
                        best_dist = self.e[v][u]
        return (best_u,best_v),best_dist

    def go_off(self):
        self.on = False

    def grow(self,cluster):
        self.v += cluster.get_v()
        for e in cluster.get_e().keys():
            if e in self.e:
                for u in cluster.get_e()[e].keys():
                    self.e[e][u] = cluster.get_e()[e][u]
            else:
                self.e[e] = cluster.get_e()[e]

    def __eq__(self, other):
        return self.name == other.name()

def heterogeneity(vertex):
    return Tree.dist[vertex].loc[vertex].sum().sum() / (len(vertex) * 2)


class Tree():
    dist = None
    n_partitions = 0
    n_cuts = 0
    he = []
    regions = []
    def __init__(self, vertex, neighbors, edges, h, middle_edge):
        self.vertex = vertex
        self.neighbors = neighbors
        self.h = h
        self.edges = {}
        if type(edges) is list:
            self.create_arch(edges)
        else:
            self.edges = edges
        self.middle_edge = middle_edge
        if self.middle_edge is None:
            middle_edge, best = None, 0
            for edge in self.edges.keys():
                middle = min(len(self.edges[edge][min(edge)]),len(self.edges[edge][max(edge)]))
                if middle > best:
                    middle_edge, best = edge, middle
            self.middle_edge = middle_edge

    def h_after_cut(self,e):
        return heterogeneity(self.edges[e][e[0]]), heterogeneity(self.edges[e][e[1]])

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

    def find_regions(self, limit_vertex, n_split, sc):
        Tree.regions += [self.vertex]
        Tree.he += [self.h]
        if len(self.vertex) > limit_vertex:
            new_trees = self.tree_partitioning(n_split,sc)
            Tree.n_partitions += 1
            for tree in new_trees:
                tree.find_regions(limit_vertex, n_split, sc)

    def tree_partitioning(self,z,sc):
        e,h,h1,h2 = self.find_best_cut(z,sc)
        new_trees = []
        Tree.n_cuts += z
        for i in range(z):
            tree_edges = self.edges.copy()
            neighbors = self.neighbors.copy()
            vertex1, vertex2 = tree_edges[e[i]][e[i][0]], tree_edges[e[i]][e[i][1]]
            edges1, edges2 = {}, {}
            neighbors[e[i][0]] = [v for v in neighbors[e[i][0]] if v!=e[i][1]]
            neighbors[e[i][1]] = [v for v in neighbors[e[i][1]] if v!=e[i][0]]
            del tree_edges[e[i]]
            middle_edge1, middle_edge2 = None, None
            best_middle1, best_middle2 = 0, 0
            for edge in tree_edges.keys():
                if edge[0] in vertex1:
                    edges1[edge] = {edge[0]:[v for v in tree_edges[edge][edge[0]] if v in vertex1],
                                    edge[1]:[v for v in tree_edges[edge][edge[1]] if v in vertex1]}
                    middle = min(len(edges1[edge][edge[0]]),len(edges1[edge][edge[1]]))
                    if middle > best_middle1:
                        best_middle1 = middle
                        middle_edge1 = edge
                else:
                    edges2[edge] = {edge[0]: [v for v in tree_edges[edge][edge[0]] if v in vertex2],
                                    edge[1]: [v for v in tree_edges[edge][edge[1]] if v in vertex2]}
                    middle = min(len(edges2[edge][edge[0]]), len(edges2[edge][edge[1]]))
                    if middle > best_middle2:
                        best_middle2 = middle
                        middle_edge2 = edge
            new_trees += [Tree(vertex1,neighbors,edges1,h1[i],middle_edge1)]
            new_trees += [Tree(vertex2, neighbors, edges2, h2[i], middle_edge2)]
        return new_trees

    def find_best_cut(self, z, sc):
        best_h1, best_h2 = [0.0]*z, [0.0]*z
        best_h1[0], best_h2[0] = self.h_after_cut(self.middle_edge)
        best_cut = [self.middle_edge]*z
        best_h = [float('inf')]*z
        best_h[0] = best_h1[0] + best_h2[0]
        max_best_h = max(best_h)
        index_max_h = best_h.index(max_best_h)
        calculated = [self.middle_edge]
        i_since_upgrade = 0
        L = {}
        to_expand = self.middle_edge
        while i_since_upgrade < sc and len(calculated) < len(self.edges):
            for v1 in to_expand:
                for v2 in self.neighbors[v1]:
                    e = (min(v1,v2),max(v1,v2))
                    if e not in calculated:
                        calculated += [e]
                        h1, h2 = self.h_after_cut(e)
                        h = h1 + h2
                        L[e] = max(h1,h2)
                        if h < max_best_h:
                            best_h[index_max_h], best_cut[index_max_h] = h, e
                            best_h1[index_max_h], best_h2[index_max_h] = h1, h2
                            max_best_h = max(best_h)
                            index_max_h = best_h.index(max_best_h)
                            i_since_upgrade = 0
            to_expand = min(L,key=L.get)
            del L[to_expand]
            i_since_upgrade += 1
        return best_cut, best_h, best_h1, best_h2

    def regions_to_onehot(self):
        R = pd.DataFrame(0,index=self.vertex,columns=range(len(Tree.regions)))
        for i in range(len(Tree.regions)):
            for el in Tree.regions[i]:
                R.loc[el,i] = 1
        return R

    def clear(self):
        Tree.regions = []
        Tree.n_cuts = 0
        Tree.n_partitions = 0
        Tree.he = []

