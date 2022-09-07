from data import *
from Solver import Solver
from LNS import LNS
from redcap import *

k=10
print('start data preparation')
t1 = t.time()
CA_vote_2004 = percentage_bush_2004()
t2 = t.time()
print('end data preparation '+ str(t2-t1))
print('number of nodes: '+ str(len(CA_vote_2004)))
root = Solver(CA_vote_2004,['bush_votes_perc'],talk=True)
res,x = root.best_k_regions(k,0.01,5,[(0,2)],[(0,100)],5)

res_redcap, x_redcap, t_redcap = redcap(CA_vote_2004,k,['bush_votes_perc'],talk=True)

print('my model total heterogeneity : ', res)

print('start build-in LNS')
t1 = t.time()
h1,r1 = root.LNS(k,10,50)
t2 = t.time()
print('end build-in LNS '+ str(t2-t1))
print('build-in LNS heterogeneity : ', sum(h1.values()))

print('start LNS')
t1 = t.time()
h2,r2 = LNS(root,k,10,50)
t2 = t.time()
print('end LNS '+ str(t2-t1))
print('LNS heterogeneity : ', sum(h2.values()))

# indicators = ['density','gdp_inhabitant','median_age','rate_migration']
# k=10
# print('start data preparation')
# t1 = t.time()
# ecodemo_eu = ecodemoEurope(indicators,2019,3)
# t2 = t.time()
# print('end data preparation '+ str(t2-t1))
# print('number of nodes: '+ str(len(ecodemo_eu)))
# root = Solver(ecodemo_eu,indicators,talk=True)
# res,x = root.best_k_regions(k,0.01,6,[(0,2)],[(0,50)],5)
# print('my model total heterogeneity : ', res)
#
#
# res_redcap, x_redcap, t_redcap = redcap(ecodemo_eu,k,indicators,talk=True)
#
# h,r=root.LNS(10,10,10)
