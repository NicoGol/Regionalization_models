import pygeoda as pg
from HierPartTree import heterogeneity
import time as t

def redcap(data,k,cols,method="fullorder-completelinkage",talk=False):
    if talk:
        print('start redcap')
    t1 = t.time()
    votes_2004 = pg.open(data)
    w = pg.queen_weights(votes_2004)
    res_redcap = pg.redcap(k,w,data[cols],method)
    t2 = t.time()
    if talk:
        print('end redcap '+ str(t2-t1))

    clusters_r = {}
    for i, c in enumerate(res_redcap['Clusters']):
        if c in clusters_r:
            clusters_r[c] += [data.index[i]]
        else:
            clusters_r[c] = [data.index[i]]
    h_r = [heterogeneity(clusters_r[i]) for i in clusters_r.keys()]
    res = sum(h_r)
    if talk:
        print('redcap total heterogeneity : ', res)
    return res, clusters_r, t2-t1
