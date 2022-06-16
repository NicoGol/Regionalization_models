import numpy as np
import pandas as pd
import time as t
from HierPartTree import HPTree, heterogeneity
from Cluster import Cluster
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


    def find_SCT_full_order_CL(self):
        i = 0
        T = {}
        edges = {}
        if self.talk:
            print('start SCT')
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

    def find_SCT_full_order_AL(self):
        i = 0
        T = {}
        edges = {}
        if self.talk:
            print('start SCT')
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
        # t1 = t.time()
        # Rs = self.MST.regions_to_onehot()
        # with gp.Env(empty=True) as env:
        #     env.setParam('OutputFlag', 0)
        #     env.start()
        #     best_obj = 0
        #     best_regions = []
        #     for i in range(len(Rs)):
        #         R = Rs[i]
        #         costs = np.array(self.MST.he[i])
        #         with gp.Model(env=env) as model:
        #             n_regions = len(self.MST.regions[i])
        #             n_points = len(self.MST.vertex)
        #             x = model.addMVar(shape=n_regions, vtype=GRB.BINARY)
        #             model.setObjective(costs @ x)
        #             model.addConstr(x.sum() == k)
        #             sparse_R = sp.csr_matrix(R.values)
        #             val = np.array([1 for _ in range(n_points)])
        #             model.addConstr(sparse_R @ x == val)
        #             model.optimize()
        #             t2 = t.time()
        #             self.times['select regions'] = t2-t1
        #             if self.talk:
        #                 print('end selecting regions ' + str(self.times['select regions']))
        #             regions = [self.MST.regions[i][j] for j,val in enumerate(x.X) if val > 0.0]
        #             if best_obj == 0 or model.objVal < best_obj:
        #                 best_obj, best_regions = model.objVal, regions
        #     return best_obj, best_regions

    def best_k_regions(self,k,limit_vertex, limit_depth, z, sc, n_best, min_size=0.0, split_tree=False):
        T = self.find_SCT_full_order_CL()
        if self.talk:
            print('start finding regions')
        t1 = t.time()
        self.SCT.find_regions(limit_vertex, limit_depth, z, sc, 0, n_best,split_tree=split_tree,min_size=min_size)
        t2 = t.time()
        self.times['column generation'] = t2-t1
        if self.talk:
            print('end finding regions ' + str(self.times['column generation']))
        res = self.select_k_regions(k)
        return res
