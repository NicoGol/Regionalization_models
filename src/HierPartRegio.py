import numpy as np
import pandas as pd
import time as t
import random
from HierPartTree import HPTree, heterogeneity, clearHPT
from Cluster import Cluster
from scipy.sparse.csgraph import minimum_spanning_tree
from Solution import Solution


class HPRegio:
    def __init__(self, contiguity_matrix, distance_matrix, talk=False):
        self.cont_m = contiguity_matrix
        self.dist = distance_matrix
        self.sols = {}
        self.vertex = self.cont_m.index.to_list()
        self.talk = talk
        self.times = {}
        self.FullO_E = None
        self.clusters = None
        clearHPT()
        HPTree.dist = self.dist
        self.dist_C = self.dist.copy()
        self.cont_C = self.cont_m.copy()
        self.SCT = None
        self.SCT_type = None

    def create_clusters(self):
        clusters = {}
        for i in self.cont_m.index:
            e = {}
            for j in self.cont_m.index:
                if i != j and self.cont_m.loc[i,j] > 0:
                    e[j] = {i: self.dist.loc[i,j]}
            clusters[i] = Cluster(i, e)
        return clusters

    def full_order_edges(self):
        Full_order_edges = {}
        indexes = self.dist.index
        for i in range(len(indexes) - 1):
            index1 = indexes[i]
            for j in range(i + 1, len(indexes)):
                index2 = indexes[j]
                dist = self.dist.loc[index1,index2]
                if index1 < index2:
                    Full_order_edges[(index1, index2)] = dist
                else:
                    Full_order_edges[(index2, index1)] = dist
        return Full_order_edges

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

    def find_MST(self):
        if self.talk:
            print('start MST')
        t1 = t.time()
        FirstO = self.cont_m.values * self.dist.values
        MST = minimum_spanning_tree(FirstO).toarray()
        neighbors = {v: [] for v in self.vertex}
        edges = {}
        for i, row in enumerate(MST):
            i = self.vertex[i]
            for j, el in enumerate(row):
                j = self.vertex[j]
                if el > 0:
                    neighbors[i].append(j)
                    neighbors[j].append(i)
                    edges[min(i,j),max(i,j)] = el
        self.SCT = HPTree(self.vertex,neighbors,list(edges.keys()),0)
        self.SCT.h = heterogeneity(self.vertex)
        self.SCT_type = 'MST'
        t2 = t.time()
        self.times['CST'] = t2 - t1
        self.times['FuO edges'] = 0.0
        self.times['Sort FuO edges'] = 0.0
        self.times['create clusters'] = 0.0
        if self.talk:
            print('end MST ' + str(self.times['CST']))
        return edges

    def find_SCT_full_order_CL(self):
        if self.FullO_E is None:
            if self.talk:
                print('start full_order')
            t1 = t.time()
            self.FullO_E = self.full_order_edges()
            t2 = t.time()
            self.times['FuO edges'] = t2-t1
            if self.talk:
                print('end full_order ' + str(self.times['FuO edges']))
                print('start sort full order')
            t1 = t.time()
            self.sort_FullO_E()
            t2 = t.time()
            self.times['Sort FuO edges'] = t2 - t1
            if self.talk:
                print('stop sort full order ' + str(self.times['Sort FuO edges']))

        if self.clusters is None:
            if self.talk:
                print('start create clusters')
            t1 = t.time()
            self.clusters = self.create_clusters()
            t2 = t.time()
            self.times['create clusters'] = t2 - t1
            if self.talk:
                print('end create clusters ' + str(self.times['create clusters']))

        i = 0
        T = {}
        edges = {v : [] for v in self.vertex}
        if self.talk:
            print('start SCT')
        t1 = t.time()
        for (u,v) in self.FullO_E.keys():
            l, m = self.clusters[u].get_name(), self.clusters[v].get_name()
            if l != m and self.cont_C[l][m] != 0 and self.FullO_E[(u,v)] >= self.dist_C[l][m]:
                e,cost = self.clusters[u].find_shortest_edge(self.clusters[v])
                T[e] = cost
                self.merge(l,m) # changing clusters, changing contiguity and changing dist_C
                edges[e[0]] += [e[1]]
                edges[e[1]] += [e[0]]
                i += 1
        t2 = t.time()
        self.times['SCT'] = t2-t1
        if self.talk:
            print('stop SCT ' + str(self.times['SCT']))
            print('start Tree arch')
        t1 = t.time()
        self.SCT = HPTree(self.vertex,edges,list(T.keys()),0,None)
        self.SCT.h = heterogeneity(self.vertex)
        self.SCT_type = 'FuO_CL'
        t2 = t.time()
        self.times['Tree arch'] = t2 - t1
        if self.talk:
            print('stop Tree arch ' + str(self.times['Tree arch']))
        return T

    def best_k_regions(self,k, limit_depth, w, sc, n_best, method='FuO_CL'):
        if method == 'FuO_CL' and self.SCT_type != method:
            T = self.find_SCT_full_order_CL()
        elif method == 'MST' and self.SCT_type != method:
            T = self.find_MST()
        if self.talk:
            print('start finding regions')
        t1 = t.time()
        self.SCT.find_regions(k, limit_depth, w, sc, n_best)
        t2 = t.time()
        self.times['regions generation'] = t2-t1
        if self.talk:
            print('end finding regions ' + str(self.times['regions generation']))
        for i in range(1,k+1):
            self.sols[i] = Solution(i,self.SCT.best_h[i][0],self.SCT.best_h[i][1])
        return self.sols[k].H, self.sols[k].to_list()

    def LNS(self,k,sc,max_iter=100):
        if self.talk:
            print('start LNS')
        t1 = t.time()
        best_sol = self.sols[k]
        current_sol = best_sol
        to_test = [(i,j) for i in range(k-1) for j in range(i+1,k)]
        sols = {}
        ok = True
        iter = 0
        visited = {}
        n_sols = 0
        while ok and iter < max_iter:
            new_sols = {}
            for (i,j) in to_test:
                r1, r2 = None, None
                ri, rj = current_sol.get_region(i), current_sol.get_region(j)
                if str(ri+rj) in visited:
                    if visited[str(ri+rj)] is not None:
                        (r1, r2) = visited[str(ri+rj)]
                elif sum(sum(self.cont_m.loc[ri,rj].values))>0:
                    r1, r2 = self.regroup_and_cut(current_sol.regions[i],current_sol.regions[j],sc)
                    visited[str(ri+ rj)] = (r1, r2)
                else:
                    visited[str(ri+ rj)] = None
                if r1 is not None and r1.h + r2.h < current_sol.get_h(i) + current_sol.get_h(j):
                    new_sol = current_sol.modify(i,j,r1,r2)
                    new_sols[(i,j)] = new_sol #sols[len(sols)] = new_sol
            if len(new_sols) > 0:
                ((i,j), min_sol) = min(new_sols.items(), key=lambda x: x[1].H)
                for (m,n) in new_sols.keys():
                    if m == i or n == i or m == j or n == j:
                        sols[n_sols] = new_sols[(m,n)]
                        n_sols += 1
            if len(sols) > 0:
                (ind,min_sol) = min(sols.items(), key= lambda x: x[1].H)
                iter += 1
                if min_sol.H < best_sol.H:
                    best_sol = min_sol
                    iter = 0
                current_sol = min_sol
                del sols[ind]
            else:
                if self.talk:
                    print('finish LNS')
                ok = False
        t2 = t.time()
        self.times['LNS'] = t2 - t1
        if self.talk:
            print('end LNS ' + str(self.times['LNS']))
        self.sols[k] = best_sol
        return self.sols[k].H, self.sols[k].to_list()

    def regroup_and_cut(self,r1,r2,sc):
        min_dist = float('inf')
        min_edge = None
        for v1 in r1.vertex:
            for v2 in r2.vertex:
                if self.cont_m.loc[v1,v2]>0 and self.dist.loc[v1,v2] < min_dist:
                    min_dist = self.dist.loc[v1,v2]
                    min_edge = (v1,v2)
        if min_edge is not None:
            (v1,v2) = min_edge
            vertex = r1.vertex + r2.vertex
            center = (min(min_edge),max(min_edge))
            edges = {}
            neighbors = {}
            for v in r1.vertex:
                neighbors[v] = r1.neighbors[v].copy()
            for v in r2.vertex:
                neighbors[v] = r2.neighbors[v].copy()
            neighbors[v1] += [v2]
            neighbors[v2] += [v1]
            for e in r1.edges.keys():
                if v1 in r1.edges[e][e[0]]:
                    edges[e] = {e[0]: r1.edges[e][e[0]].copy() + r2.vertex, e[1]: r1.edges[e][e[1]].copy()}
                else:
                    edges[e] = {e[0]: r1.edges[e][e[0]].copy(), e[1]: r1.edges[e][e[1]].copy() + r2.vertex}
            for e in r2.edges.keys():
                if v2 in r2.edges[e][e[0]]:
                    edges[e] = {e[0]: r2.edges[e][e[0]].copy() + r1.vertex, e[1]: r2.edges[e][e[1]].copy()}
                else:
                    edges[e] = {e[0]: r2.edges[e][e[0]].copy(), e[1]: r2.edges[e][e[1]].copy() + r1.vertex}
            if center[0] in r1.vertex:
                edges[center] = {center[0]:r1.vertex, center[1]:r2.vertex}
            else:
                edges[center] = {center[0]: r2.vertex, center[1]: r1.vertex}
            concatree = HPTree(vertex,neighbors,edges,0, center={center:(0,0)})
            [new_r1, new_r2] = concatree.tree_partitioning(1,sc,1)
            if new_r1.vertex == r1.vertex or new_r1.vertex == r2.vertex:
                new_r1, new_r2 = None, None
            return new_r1, new_r2

