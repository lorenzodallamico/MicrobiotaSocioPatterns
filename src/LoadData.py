import pandas as pd
import networkx as nx
import numpy as np
from skbio import TreeNode


class ReturnValue:
    def __init__(self, df_net, G, microbiota_all, asv_ids, tree, days):
        self.df_net = df_net
        self.G = G
        self.microbiota_all = microbiota_all
        self.asv_ids = asv_ids
        self.tree = tree
        self.days = days

def Load(filename):
    '''Loads and assembles the proximity network and the microbiome dataset used throughout the
    analysis: the contact network (from the two SocioPatterns deployments), the ASV taxonomy table,
    the rooted phylogenetic tree needed for UniFrac, and the microbiota abundance table restricted
    to the first 5 sampling days.

    Use: data = Load(filename)

    Inputs:
        * filename (str): path to the cleaned microbiota csv (e.g. 'Data/ASV_cleaned.csv')
    '''
    
    df_net = pd.read_csv('Data/df_net.csv', dtype={'pid': str, 'pid2': str})
    df_net = df_net.groupby(['pid', 'pid2']).weight.sum().reset_index()

    G = nx.from_pandas_edgelist(df_net, source = "pid", target = "pid2", edge_attr = ['weight'])

    asv_ids = pd.read_excel('Data/table__tax_by_ASV.xlsx')['Taxonomy of each ASV'].iloc[4:].values

    # UniFrac requires a rooted tree
    tree = TreeNode.read('Data/molde__nwk_trees.tar-1/6-rooted-tree_nwk/6-rooted-tree.nwk')

    # microbiome data, restricted to the first 5 sampling days
    df_microbiota = pd.read_csv(f'{filename}', dtype={'Person ID': str})

    d_th = np.sort(df_microbiota.Day.unique())[4]
    df_microbiota = df_microbiota[df_microbiota.Day <= d_th]

    meta_cols = [c for c in ['Sample ID', 'Sampling-date', 'Day'] if c in df_microbiota.columns]
    microbiota_all = df_microbiota.set_index(['Person ID', 'Day code']).drop(columns=meta_cols)

    days = np.sort(np.unique([x[1] for x in microbiota_all.index]))

    return ReturnValue(df_net, G, microbiota_all, asv_ids, tree, days)
