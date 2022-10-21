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
