import pandas as pd
import numpy as np

def estimate_richness(df, groupby='bindCol', query="fragmentCountAln_adj > 0"):
    richness = df.query(query).groupby(groupby).agg(
        {
            'gene': 'nunique', 
            'fragmentCountAln_adj': 'sum'
        }
    ).rename(columns={'gene': 'n_references'})


    zero_richness = pd.DataFrame(
        {
            groupby: df.loc[~df[groupby].isin(richness.index), groupby].unique(),
            'n_references': 0,
            'fragmentCountAln_adj': 0
        },
    ).dropna().set_index(groupby)
    #zero_richness.index.name = groupby

    return  pd.concat([richness, zero_richness])

def estimate_chao1(df, groupby=['bindCol','gene'], query="fragmentCountAln_adj > 0"):
    gene_counts = df.query(query).groupby(groupby).agg(
        {
            'fragmentCountAln_adj': 'sum', 
        }
    ).reset_index()

    N = gene_counts.groupby(groupby[0]).agg({groupby[1]: 'nunique'}).rename(columns={groupby[1]: 'N'})
    S = gene_counts.query("fragmentCountAln_adj == 1").groupby(groupby[0]).agg({groupby[1]: 'nunique'}).rename(columns={groupby[1]: 'S'})
    D = gene_counts.query("fragmentCountAln_adj == 2").groupby(groupby[0]).agg({groupby[1]: 'nunique'}).rename(columns={groupby[1]: 'D'})

    chao1 = pd.concat([N, S, D], axis=1).fillna(0)
    chao1['chao1'] = chao1.apply(lambda x: x['N'] + x['S']**2 / (2 * x['D']), axis=1)
    chao1['chao1_bc'] = chao1.apply(lambda x: x['N'] + x['S']*(x['S']-1) / (2 * (x['D']+1)), axis=1)

    return chao1


def estimate_shannon(df, groupby=['bindCol', 'gene'], query="fragmentCountAln_adj > 0", totcol=None):
    shannon = df.query(query).groupby(groupby).agg(
        {
            'fragmentCountAln_adj': 'sum'
        }
    ).reset_index()
    
    if totcol is None:
        totcol = df.query(query).groupby(groupby[0]).agg({ 'fragmentCountAln_adj': 'sum'}).reset_index()
        
    shannon = shannon.merge(
        totcol,
        on=groupby[0],
        suffixes=('_gene', '_total')
    )
    shannon['p'] = shannon['fragmentCountAln_adj_gene'] / shannon['fragmentCountAln_adj_total']
    shannon['eq'] = shannon['p'] * np.log(shannon['p'])
    shannon = -shannon.groupby(groupby[0]).agg({'eq': 'sum'}).rename(columns={'eq': 'shannon'})
        

    return  shannon

def estimate_simpson(df, groupby=['bindCol', 'gene'], query="fragmentCountAln_adj > 0", totcol=None):
    simpson = df.query(query).groupby(groupby).agg(
        {
            'fragmentCountAln_adj': 'sum'
        }
    ).reset_index()
    
    if totcol is None:
        totcol = df.query(query).groupby(groupby[0]).agg({ 'fragmentCountAln_adj': 'sum'}).reset_index()
        
    simpson = simpson.merge(
        totcol,
        on=groupby[0],
        suffixes=('_gene', '_total')
    )
    simpson['p'] = simpson['fragmentCountAln_adj_gene'] / simpson['fragmentCountAln_adj_total']
    simpson['eq'] = simpson['p'] ** 2
    simpson = simpson.groupby(groupby[0]).agg({'eq': 'sum'}).rename(columns={'eq': 'simpson'})
        

    return  simpson