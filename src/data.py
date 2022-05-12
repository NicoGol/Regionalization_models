import pandas as pd
import geopandas


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


