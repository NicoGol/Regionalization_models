import pandas as pd
import geopandas
from sklearn.preprocessing import MinMaxScaler


def percentage_bush_2004(states='all'):
    df = pd.read_csv('../data/countypres_2000-2020.csv')
    df_2004 = df.loc[(df.year==2004) & (df.party=='REPUBLICAN')]
    if states != 'all':
        df_2004 = df_2004.loc[df_2004.state.isin(states)]
    df_2004['bush_votes_perc'] = df_2004['candidatevotes'] / df_2004['totalvotes']
    df_2004 = df_2004.drop(df_2004[df_2004.bush_votes_perc.isnull()].index)
    df_2004['county_fips'] = df_2004['county_fips'].apply(int)

    contiguity = geopandas.read_file('../data/ncovr/NAT.shp')
    contiguity['FIPS'] = contiguity['FIPS'].apply(int)

    gdf = contiguity.merge(df_2004, left_on='FIPS', right_on='county_fips')
    gdf['neighbors'] = [[]]*len(gdf)
    for index, county in gdf.iterrows():
        # get 'not disjoint' countries
        neighbors = gdf[~gdf.geometry.disjoint(county.geometry)].county_fips.tolist()
        neighbors.remove(county.county_fips)

        # add names of neighbors as NEIGHBORS value
        gdf.at[index, 'neighbors'] = neighbors
    #gdf = gdf.set_index('county_fips',drop=True)#[['county_fips', 'county_name', 'state', 'state_po', 'bush_votes_perc','candidatevotes','totalvotes','neighbors','geometry']].set_index('county_fips',drop=True)
    return gdf

def ecodemoEurope(indicators,year,level):
    df = None
    for ind in indicators:
        df_ind = pd.read_csv('../data/EcoDemoEurope/'+ind+'.csv')[['geo','TIME_PERIOD','OBS_VALUE']]
        df_ind.rename(columns={'OBS_VALUE':ind},inplace=True)
        if df is None:
            df = df_ind
        else:
            df = pd.merge(left=df,right=df_ind,on=['geo','TIME_PERIOD'],how='inner')
    df = df[df['TIME_PERIOD']==year]
    df = df[df['geo'].apply(lambda row: row[:2] != 'CY' and row[:2] != 'IE' and row[:2] != 'MT')]
    #df = df[df['geo'].apply(lambda row: len(row)==2+block_size)]
    gdf = geopandas.read_file('../data/EcoDemoEurope/NUTS_RG_20M_2021_3035.shp')
    gdf = gdf[gdf['LEVL_CODE']==level]
    gdf = gdf.merge(df,left_on='NUTS_ID',right_on='geo')
    gdf['neighbors'] = [[] for _ in range(len(gdf))]
    for index, row in gdf.iterrows():
        # get 'not disjoint' countries
        neighbors = gdf[~gdf.geometry.disjoint(row.geometry)].NUTS_ID.tolist()
        neighbors.remove(row.NUTS_ID)

        # add names of neighbors as NEIGHBORS value
        gdf.at[index, 'neighbors'] = neighbors
    scaler = MinMaxScaler()
    gdf[indicators] = scaler.fit_transform(gdf[indicators])
    gdf = gdf.set_index('NUTS_ID',drop=True)
    visited = []
    to_expand = [gdf.index[0]]
    while to_expand:
        next = to_expand[0]
        if next not in visited:
            visited.append(next)
            to_expand += gdf['neighbors'][next]
        to_expand = to_expand[1:]
    gdf = gdf.loc[visited]
    return gdf



