from data import *
from HPR import HPR
from literature_models import *
import cls

# k=5
# edu_BE = pd.read_pickle('../data/education_BE.pkl')
# edu_BE_dist = pd.read_pickle('../data/education_BE_dist.pkl')
# edu_BE_cont = pd.read_pickle('../data/education_BE_cont.pkl')
# HPR = HPR(k,[(0,2)],edu_BE_cont,edu_BE_dist,LNS=False,talk=True)
# H_redcap, regions_redcap, t_redcap = redcap(edu_BE,k,['EDU_LOW_r', 'EDU_MID_r', 'EDU_HIGH_r'],edu_BE_dist)

############################################################################

# k=15
# Ecoregions_dist = pd.read_pickle('../data/USA_ecoregions_dist.pkl')
# Ecoregions_cont = pd.read_pickle('../data/USA_ecoregions_cont.pkl')
# HPR,h,regions = HPR(k,[(0,2)],Ecoregions_cont,Ecoregions_dist,LNS=True,talk=True)

############################################################################

# EU_NUTS1 = pd.read_pickle('../data/ecodemo_NUTS1.pkl')
# EU_NUTS1_cont = pd.read_pickle('../data/ecodemo_NUTS1_cont.pkl')
# EU_NUTS1_dist = pd.read_pickle('../data/ecodemo_NUTS1_dist.pkl')
# k=5
# HPR = HPR(k,[(0,2)],EU_NUTS1_cont,EU_NUTS1_dist,LNS=False,talk=True)
# H_redcap, regions_redcap, t_redcap = redcap(EU_NUTS1,k,['density','gdp_inhabitant','median_age','rate_migration'],EU_NUTS1_dist)

#############################################################################

# EU_NUTS2 = pd.read_pickle('../data/ecodemo_NUTS2.pkl')
# EU_NUTS2_cont = pd.read_pickle('../data/ecodemo_NUTS2_cont.pkl')
# EU_NUTS2_dist = pd.read_pickle('../data/ecodemo_NUTS2_dist.pkl')
# k=5
# HPR = HPR(k,[(0,2)],EU_NUTS2_cont,EU_NUTS2_dist,max_depth=4,LNS=False,talk=True)
# H_redcap, regions_redcap, t_redcap = redcap(EU_NUTS2,k,['density','gdp_inhabitant','median_age','rate_migration'],EU_NUTS2_dist)

#############################################################################

# EU_NUTS3_cont = pd.read_pickle('../data/NUTS3_all_cont.pkl')
# EU_NUTS3_dist = pd.read_pickle('../data/NUTS3_all_dist.pkl')
# k=10
# HPR,h,regions = HPR(k,[(0,2)],EU_NUTS3_cont,EU_NUTS3_dist,LNS=True,talk=True)

############################################################################

# k=10
# USA_vote_cont = pd.read_pickle('../data/USA_vote_2004_cont.pkl')
# USA_vote_dist = pd.read_pickle('../data/USA_vote_2004_dist.pkl')
# HPR,h,regions = HPR(k,[(0,2)],USA_vote_cont,USA_vote_dist,LNS=True,talk=True)

############################################################################

# k=10
# print('start data preparation')
# t1 = t.time()
# CA_vote_2004 = percentage_bush_2004()
# t2 = t.time()
# print('end data preparation '+ str(t2-t1))
# print('number of nodes: '+ str(len(CA_vote_2004)))
# root = Solver(CA_vote_2004,['bush_votes_perc'],talk=True)
# res,x = root.best_k_regions(k,0.01,5,[(0,2)],[(0,100)],5)
#
# res_redcap, x_redcap, t_redcap = redcap(CA_vote_2004,k,['bush_votes_perc'],talk=True)
#
# print('my model total heterogeneity : ', res)
#
# print('start build-in LNS')
# t1 = t.time()
# h1,r1 = root.LNS(k,10,50)
# t2 = t.time()
# print('end build-in LNS '+ str(t2-t1))
# print('build-in LNS heterogeneity : ', sum(h1.values()))
#
# print('start LNS')
# t1 = t.time()
# h2,r2 = LNS(root,k,10,50)
# t2 = t.time()
# print('end LNS '+ str(t2-t1))
# print('LNS heterogeneity : ', sum(h2.values()))


##############################################################################################

#indicators = ['density','gdp_inhabitant','median_age','rate_migration']
# k=10
# print('start data preparation')
# t1 = t.time()
#ecodemo_eu = ecodemoEurope(indicators,2019,3)
#dist = generate_dist_matrix(ecodemo_eu,indicators)
#cont = generate_cont_matrix(ecodemo_eu)
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


