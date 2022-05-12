import gurobipy as gp
from gurobipy import GRB
from data import *
from Cluster import *
import time as t
from redcap import *

k=10
print('start data preparation')
t1 = t.time()
CA_vote_2004 = percentage_bush_2004()#['CALIFORNIA','NEVADA'])
t2 = t.time()
print('end data preparation '+ str(t2-t1))
print('number of nodes: '+ str(len(CA_vote_2004)))
#root = Graph(CA_vote_2004.set_index('county_fips',drop=True),['bush_votes_perc'],talk=True)
#res,x = root.best_k_regions(k,10,[(0,2)],10,split_tree=True)

res_redcap, t_redcap = redcap(CA_vote_2004,k,talk=True)

#print('my model total heterogeneity : ', res)
