from data import *
import matplotlib.pyplot as plt
import pygeoda
import time as t


df = percentage_bush_2004(state='CALIFORNIA')
print(len(df))
x = range(len(df))
perc_votes = list(df.sort_values(by=['bush_votes_perc']).bush_votes_perc)
df_sort = df.sort_values(by=['totalvotes'])
total_votes = list(df_sort.totalvotes)
bush_votes = list(df_sort.candidatevotes)
plt.plot(x,perc_votes)
plt.title('Percentage of votes for Bush in 2004 per county')
plt.show()
plt.plot(x,total_votes,'b',x,bush_votes,'r')
plt.title('Number of votes for US 2004 election per county')
plt.yscale('log')
plt.show()


#w1 = pygeoda.queen_weights(df)
#print(w1)
#w2 = pygeoda.rook_weights(df)
#print(w2)
