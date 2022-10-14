import numpy as np
import pandas as pd
import time as t
import random
from HierPartTree import HPTree, heterogeneity
from Cluster import Cluster
import gurobipy as gp
from gurobipy import GRB
import scipy.sparse as sp

class Solver:
    def __init__(self, df, cols,vertex=None,first_order=None,full_order=None, dist_deg=2, talk=False):
        self.df = df
        self.cols = cols
        self.regions = {}
        self.he = {}
        self.H = {}
        if vertex is None:
            self.vertex = df.index.to_list()
        else:
            self.vertex = vertex
        self.talk = talk
        self.times = {}
        self.FirstO_E, self.cont_C, self.cont_m, self.clusters = None, None, None, None
        if first_order is None:
            if talk:
                print('start first_order')
            t1 = t.time()
            self.first_order_edges(df, cols, dist_deg)
            t2 = t.time()
            self.times['FiO edges'] = t2-t1
        else:
            self.FirstO_E, self.cont_C, self.cont_m, self.clusters = first_order[0], first_order[1].copy(), first_order[1], first_order[2]
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
        HPTree.dist = self.dist
        self.dist_C = self.dist.copy()
        self.SCT = None

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
        self.FirstO_E, self.cont_m, self.cont_C, self.clusters = First_order_edges, contiguity_matrix, contiguity_matrix.copy(), clusters

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
        self.cont_C[l] = self.cont_C[[l,m]].max(axis=1)
        self.cont_C.loc[l] = self.cont_C.loc[[l, m]].max(axis=0)
        self.dist_C[l] = self.dist_C[[l, m]].max(axis=1)
        self.dist_C.loc[l] = self.dist_C.loc[[l, m]].max(axis=0)


    def find_SCT_full_order_CL(self):
        i = 0
        T = {}
        edges = {}
        if self.talk:
            print('start SCT')
        t1 = t.time()
        for (u,v) in self.FullO_E.keys():
            l, m = self.clusters[u].get_name(), self.clusters[v].get_name()
            if l != m and self.cont_C[l][m] != 0 and self.FullO_E[(u,v)] >= self.dist_C[l][m]:
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
        self.times['SCT'] = t2-t1
        if self.talk:
            print('stop SCT ' + str(self.times['SCT']))
            print('start Tree arch')
        t1 = t.time()
        self.SCT = HPTree(self.vertex,edges,list(T.keys()),0,None)
        self.SCT.h = heterogeneity(self.vertex)
        t2 = t.time()
        self.times['Tree arch'] = t2 - t1
        if self.talk:
            print('stop Tree arch ' + str(self.times['Tree arch']))
        return T

    def select_k_regions(self,k):
        if self.talk:
            print('start selecting regions')
        t1 = t.time()
        R = self.SCT.regions_to_onehot()
        with gp.Env(empty=True) as env:
            env.setParam('OutputFlag', 0)
            env.start()
            costs = np.array(list(self.SCT.dict_he.values()))
            with gp.Model(env=env) as model:
                n_regions = len(self.SCT.dict_regions)
                n_points = len(self.SCT.vertex)
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
                regions = [self.SCT.dict_regions[list(self.SCT.dict_regions.keys())[j]] for j,val in enumerate(x.X) if val > 0.0]
                return model.objVal, regions

    def best_k_regions(self,k,limit_vertex, limit_depth, w, sc, n_best):
        if self.SCT is None:
            T = self.find_SCT_full_order_CL()
        if self.talk:
            print('start finding regions')
        t1 = t.time()
        self.SCT.find_regions(k,limit_vertex, limit_depth, w, sc, n_best,talk=self.talk)
        t2 = t.time()
        self.times['column generation'] = t2-t1
        if self.talk:
            print('end finding regions ' + str(self.times['column generation']))
        h = self.SCT.sub_h[k][0]
        regions = [HPTree.dict_regions[i] for i in self.SCT.sub_h[k][1]]
        self.regions[k] = regions
        self.he[k] = [HPTree.dict_he[i] for i in self.SCT.sub_h[k][1]]
        self.H[k] = h
        return h, regions
        # res = self.select_k_regions(k)
        # return res

    def LNS(self,k,limit,sc):
        best_h = {i: HPTree.dict_he[r] for i, r in enumerate(self.SCT.sub_h[k][1])}
        best_regions = {i:HPTree.dict_regions[r] for i,r in enumerate(self.SCT.sub_h[k][1])}
        best_neighbors = {i: {v:[e for e in HPTree.dict_neighbors[r][v] if e in best_regions[i]] for v in HPTree.dict_neighbors[r].keys() if v in best_regions[i]} for i, r in enumerate(self.SCT.sub_h[k][1])}
        best_edges = {i:[e for e in self.SCT.edges.keys() if e[0] in best_regions[i] and e[1] in best_regions[i]] for i in range(k)}
        tested = []
        not_contiguous = []
        i_since_last_up = 0
        n_comb = sum([l for l in range(k)])
        while i_since_last_up < limit and len(tested)+len(not_contiguous) < n_comb:
            l = random.sample(range(k), 2)
            l.sort()
            [i,j] = l
            if (i,j) not in tested and (i,j) not in not_contiguous:
                if sum(sum(self.cont_m.loc[best_regions[i],best_regions[j]].values))>0:
                    e, h, h1, h2, vertex1, vertex2, neighbors, edges = self.regroup_and_cut(best_regions[i],best_regions[j],best_neighbors[i],best_neighbors[j],best_edges[i],best_edges[j],sc)
                    if h < best_h[i] + best_h[j]:
                        neighbors[e[0]] = [v for v in neighbors[e[0]] if v != e[1]]
                        neighbors[e[1]] = [v for v in neighbors[e[1]] if v != e[0]]
                        best_h[i] = h1
                        best_h[j] = h2
                        best_regions[i] = vertex1
                        best_regions[j] = vertex2
                        best_neighbors[i] = {v:neighbors[v] for v in neighbors.keys() if v in vertex1}
                        best_neighbors[j] = {v:neighbors[v] for v in neighbors.keys() if v in vertex2}
                        best_edges[i] = [e for e in edges if e[0] in vertex1 and e[1] in vertex1]
                        best_edges[j] = [e for e in edges if e[0] in vertex2 and e[1] in vertex2]
                        i_since_last_up = 0
                        tested = [e for e in tested if i not in e and j not in e]
                        not_contiguous = [e for e in not_contiguous if i not in e and j not in e]
                    else:
                        i_since_last_up += 1
                    tested += [(i, j)]
                else:
                    not_contiguous += [(i,j)]
        self.regions[k] = best_regions.values()
        self.he[k] = best_h.values()
        self.H[k] = sum(self.he)
        return best_h, best_regions

    def regroup_and_cut(self,r1,r2,n1,n2,e1,e2,sc):
        min_dist = float('inf')
        min_edge = None
        for v1 in r1:
            for v2 in r2:
                if self.cont_m.loc[v1,v2]>0 and self.dist.loc[v1,v2] < min_dist:
                    min_dist = self.dist.loc[v1,v2]
                    min_edge = (v1,v2)
        if min_edge is not None:
            (v1,v2) = min_edge
            vertex = r1 + r2
            center = (min(min_edge),max(min_edge))
            neighbors = {v:n1[v].copy() for v in n1.keys()}
            for v in r2:
                neighbors[v] = n2[v].copy()
            neighbors[v1] += [v2]
            neighbors[v2] += [v1]
            edges = e1.copy() + e2.copy()#[e for e in self.SCT.edges.keys() if (e[0] in r1 and e[1] in r1) or (e[0] in r2 and e[1] in r2)]
            edges += [(v1,v2)]
            concatree = HPTree(vertex,neighbors,edges,0, center={center:(0,0)},index=-1)
            e, h, h1, h2 = concatree.find_best_cut(1,sc)
            vertex1 = concatree.edges[e[0]][e[0][0]]
            vertex2 = concatree.edges[e[0]][e[0][1]]
            return e[0], h[0], h1[0], h2[0], vertex1, vertex2, neighbors, edges

    def regions_to_gdp(self,k,colors):
        if k in self.regions:
            self.df['colors'] = None
            for i,r in enumerate(self.regions[k]):
                for v in r:
                    self.df.loc[v, 'colors'] = colors[i]
            return self.df
        else:
            return None
