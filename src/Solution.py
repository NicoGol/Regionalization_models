from Solver import *

class Solution:
    def __init__(self,k,H,regions):
        self.k = k
        if type(regions) is list:
            self.regions = {i : regions[i] for i in range(k)}
        else:
            self.regions = regions
        self.H = H

    def get_region(self,i):
        return self.regions[i].vertex

    def get_neighbors(self,i):
        return self.regions[i].neighbors

    def get_edges(self,i):
        return self.regions[i].edges

    def get_h(self,i):
        return self.regions[i].h

    def copy(self):
        return Solution(self.k,self.H,{i:r.copy() for (i,r) in self.regions.items()})

    def modify(self,i,j,r1,r2):
        new_sol = self.copy()
        new_sol.regions[i], new_sol.regions[j] = r1, r2
        new_sol.H += r1.h + r2.h - (self.get_h(i) + self.get_h(j))
        return new_sol

    def to_list(self):
        R = []
        for region in self.regions.values():
            R.append(region.vertex)
        return R


# def LNS(solver,k,limit,sc):
#     best_h = {i: HPTree.dict_he[r] for i, r in enumerate(solver.SCT.sub_h[k][1])}
#     best_regions = {i: HPTree.dict_regions[r] for i, r in enumerate(solver.SCT.sub_h[k][1])}
#     tested = []
#     not_contiguous = []
#     i_since_last_up = 0
#     n_comb = sum([l for l in range(k)])
#     while i_since_last_up < limit and len(tested) + len(not_contiguous) < n_comb:
#         l = random.sample(range(k), 2)
#         l.sort()
#         [i, j] = l
#         if (i, j) not in tested and (i, j) not in not_contiguous:
#             if sum(sum(solver.cont_m.loc[best_regions[i], best_regions[j]].values)) > 0:
#                 vertex = best_regions[i] + best_regions[j]
#                 new_solver = regroup(solver,vertex)
#                 T = new_solver.find_SCT_full_order_CL()
#                 e, h, h1, h2 = new_solver.SCT.find_best_cut(1,sc,n_best=1)
#                 e, h, h1, h2 = e[0], h[0], h1[0], h2[0]
#                 if h < best_h[i] + best_h[j]:
#                     vertex1 = new_solver.SCT.edges[e][e[0]]#list(set(vertex).intersection(new_solver.SCT.edges[e][e[0]]))
#                     vertex2 = new_solver.SCT.edges[e][e[1]]#list(set(vertex).intersection(new_solver.SCT.edges[e][e[1]]))
#                     best_h[i] = h1
#                     best_h[j] = h2
#                     best_regions[i] = vertex1
#                     best_regions[j] = vertex2
#                     i_since_last_up = 0
#                     tested = [e for e in tested if i not in e and j not in e]
#                     not_contiguous = [e for e in not_contiguous if i not in e and j not in e]
#                 else:
#                     i_since_last_up += 1
#                 tested += [(i, j)]
#             else:
#                 not_contiguous += [(i, j)]
#     return best_h, best_regions
#
# def regroup(solver,vertex):
#     FullO_E = {e:solver.FullO_E[e] for e in solver.FullO_E.keys() if e[0] in vertex and e[1] in vertex}
#     dist = solver.dist.loc[vertex,vertex]
#     df = solver.df.loc[vertex]
#     clusters = {}
#     FirstO_E = {}
#     for e in solver.FirstO_E.keys():
#         if e[0] in vertex and e[1] in vertex:
#             FirstO_E[e] = solver.FirstO_E[e]
#             if e[0] in clusters:
#                 clusters[e[0]].e[e[1]] = {e[0]:solver.FirstO_E[e]}
#             else:
#                 clusters[e[0]] = Cluster(e[0],{e[1]:{e[0]:solver.FirstO_E[e]}})
#             if e[1] in clusters:
#                 clusters[e[1]].e[e[0]] = {e[1]:solver.FirstO_E[e]}
#             else:
#                 clusters[e[1]] = Cluster(e[1],{e[0]:{e[1]:solver.FirstO_E[e]}})
#     cont_m = solver.cont_m.loc[vertex,vertex]
#     return Solver(df,solver.cols,vertex=vertex,first_order=(FirstO_E,cont_m,clusters),full_order=(FullO_E,dist))
