import pygeoda as pg
from Cluster import *

def redcap(data,k,method="fullorder-completelinkage",talk=False):
    if talk:
        print('start redcap')
    t1 = t.time()
    votes_2004 = pg.open(data)
    w = pg.queen_weights(votes_2004)
    res_redcap = pg.redcap(k,w,data[['bush_votes_perc']],method)
    t2 = t.time()
    if talk:
        print('end redcap '+ str(t2-t1))

    clusters_r = {}
    for i, c in enumerate(res_redcap['Clusters']):
        if c in clusters_r:
            clusters_r[c] += [data.loc[i,'FIPS']]
        else:
            clusters_r[c] = [data.loc[i, 'FIPS']]
    h_r = [heterogeneity(clusters_r[i]) for i in clusters_r.keys()]
    res = sum(h_r)
    if talk:
        print('redcap total heterogeneity : ', res)
    return res, clusters_r, t2-t1
