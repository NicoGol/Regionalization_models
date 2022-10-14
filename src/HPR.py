from HierPartRegio import *
import math

def HPR(k,w,contiguity_matrix,distance_matrix,LNS=False,sc=None,max_depth=None,n_best=5,method='FuO_CL',LNS_iter=20,talk=False):
    #if list(contiguity_matrix.index != distance_matrix.index or contiguity_matrix.columns != distance_matrix.columns:
        #return 'wrong parameters'
    solver = HPRegio(contiguity_matrix,distance_matrix,talk=talk)
    if max_depth is None:
        md = {5:3,10:6,15:6}
        if k in md:
            max_depth = md[k]
        else:
            max_depth = int(math.log(k,2)) +2
    if sc is None:
        sc = [(0,int(len(contiguity_matrix)/10))]
    h, regions = solver.best_k_regions(k,max_depth,w,sc,n_best,method=method)
    if LNS:
        h_LNS, regions_LNS = solver.LNS(k,sc[0][1],max_iter=LNS_iter)
    return solver