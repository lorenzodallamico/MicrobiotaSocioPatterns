import umap
import numpy as np
from scipy.stats import spearmanr
from copy import copy
from scipy.stats import mannwhitneyu
from itertools import combinations
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm

from src.distances import *


def People_vs_time(df, distance):
    '''Computes pairwise microbiota distances within individuals (across time)
    and between individuals.

    Use: sim, nsim = People_vs_time(df, distance)

    Inputs:
        * df (pandas DataFrame): abundance table with a 'Person ID' column;
          taxa abundances start from column index 4
        * distance (str): distance measure to use; one of 'jaccard',
          'bray-curtis', or 'TaxonomicJaccard'

    Outputs:
        * sim (list): distances between repeated measurements of the same individual
        * nsim (list): distances between measurements of different individuals
    '''
    sim, nsim = [], []
    all_id = set(df['Person ID'])

    if distance == 'TaxonomicJaccard':
        D = getBranches(df.iloc[:, 4:])
        indices = np.array(list(D.keys()))

    for idx in all_id:
        bool_idx = df['Person ID'] == idx
        X = df[bool_idx].iloc[:, 4:].values
        m, _ = X.shape

        for i in range(m):
            for j in range(i + 1, m):
                if distance == 'jaccard':
                    sim.append(jaccard_similarity(X[i], X[j]))
                elif distance == 'TaxonomicJaccard':
                    sim.append(TaxonomicJaccard(D, indices[bool_idx][i], indices[bool_idx][j]))
                elif distance == 'bray-curtis':
                    sim.append(bray_curtis(X[i], X[j]))

        Y = df[~bool_idx].iloc[:, 4:].values
        q, _ = Y.shape

        for i in range(m):
            for j in range(q):
                if distance == 'jaccard':
                    nsim.append(jaccard_similarity(X[i], Y[j]))
                elif distance == 'TaxonomicJaccard':
                    nsim.append(TaxonomicJaccard(D, indices[bool_idx][i], indices[~bool_idx][j]))
                elif distance == 'bray-curtis':
                    nsim.append(bray_curtis(X[i], Y[j]))

    return sim, nsim


def Update_dfnet(data):
    '''Updates data.df_net with all node pairs, zero-filling missing edges,
    and appends pairwise Jaccard and TaxonomicJaccard microbiota similarity columns.

    All individuals present in both the contact network and the microbiota
    dataset are included. Pairs with no observed contact receive weight = 0.

    Inputs:
        * data (DATA): data container with attributes df_net, G,
          and df_microbiota_bool

    Outputs:
        * (modifies data.df_net in place)
    '''
    idx1, idx2, weight = [], [], []
    all_nodes = np.unique(data.df_net[['pid', 'pid2']])
    all_nodes = all_nodes[np.isin(all_nodes, data.df_microbiota_bool.index)]

    for pid, pid2 in combinations(all_nodes, 2):
        idx1.append(pid)
        idx2.append(pid2)
        if tuple([pid, pid2]) in data.G.edges:
            weight.append(data.G[pid][pid2]['weight'])
        else:
            weight.append(0)

    df_net = pd.DataFrame(np.array([idx1, idx2, weight]).T, columns=['pid', 'pid2', 'weight'])

    df_net['jaccard'] = df_net.apply(
        lambda x: jaccard_similarity(
            data.df_microbiota_bool.loc[x.pid].values,
            data.df_microbiota_bool.loc[x.pid2]
        ), axis=1
    )

    D = getBranches(data.df_microbiota_bool)
    df_net['TaxonomicJaccard'] = df_net.apply(lambda x: TaxonomicJaccard(D, x.pid, x.pid2), axis=1)

    data.df_net = df_net
    data.df_net.weight = data.df_net.weight.astype(float)


def df_net_strong_ties(data):
    '''Generates a null model by randomly shuffling edge weights and computing
    average microbiota similarity above a series of contact-duration thresholds.

    Use: Jaccard, UniSim, tv = df_net_strong_ties(data)

    Inputs:
        * data (DATA): data container whose df_net has 'weight', 'jaccard',
          and 'TaxonomicJaccard' columns

    Outputs:
        * Jaccard (list of lists): for each randomisation trial, the mean
          Jaccard similarity above each threshold in tv
        * UniSim (list of lists): same for TaxonomicJaccard similarity
        * tv (array): contact-duration thresholds in hours (6 values from 0 to 6)
    '''
    tv = np.linspace(0, 6, 6)

    df_net_rdn = copy(data.df_net)
    x = copy(df_net_rdn.weight.values)
    Jaccard, UniSim = [], []
    n_sim = 1500

    for i in range(n_sim):
        print(i, end='\r')

        np.random.shuffle(x)
        df_net_rdn.weight = x
        Jaccard.append([df_net_rdn[df_net_rdn.weight * 10 / 3600 > t].jaccard.mean() for t in tv])
        UniSim.append([df_net_rdn[df_net_rdn.weight * 10 / 3600 > t].TaxonomicJaccard.mean() for t in tv])

    return Jaccard, UniSim, tv


def GetTaxaPValues(data):
    '''For each taxon, tests whether pairs of individuals who both carry it
    have significantly different contact duration compared with pairs where
    only one individual carries it (Mann-Whitney U test).

    Use: Psmall, Plarge, Ysingle, Yboth, rsmall, rlarge, delta_t = GetTaxaPValues(data)

    Inputs:
        * data (DATA): data container with attributes G, df_net, and df_microbiota_bool

    Outputs:
        * Psmall (dict): p-values for the alternative 'less' (both < one)
        * Plarge (dict): p-values for the alternative 'greater' (both > one)
        * Ysingle (dict): contact-duration arrays for pairs where only one
          individual carries the taxon
        * Yboth (dict): contact-duration arrays for pairs where both
          individuals carry the taxon
        * rsmall, rlarge (dict): rank-biserial correlation coefficients
        * delta_t (dict): median(both) - median(single) for each taxon
    '''
    all_nodes = np.unique(data.df_net[['pid', 'pid2']])
    idx1, idx2, weights = [], [], []

    for pid, pid2 in combinations(all_nodes, 2):
        idx1.append(pid)
        idx2.append(pid2)
        if tuple([pid, pid2]) in data.G.edges:
            weights.append(data.G[pid][pid2]['weight'])
        else:
            weights.append(0)

    df_net = pd.DataFrame(np.array([idx1, idx2, weights]).T, columns=['pid', 'pid2', 'weight'])
    df_net.weight = df_net.weight.astype(float)

    A = data.df_microbiota_bool.loc[df_net.pid].values
    B = data.df_microbiota_bool.loc[df_net.pid2].values
    IDXboth = (A * B) > 0
    IDXone = (A + B) * (1 - A * B) > 0
    _, m = IDXboth.shape

    Yboth = [df_net.weight.values[IDXboth[:, i]] for i in range(m)]
    Ysingle = [df_net.weight.values[IDXone[:, i]] for i in range(m)]

    index, Psmall, Plarge, rsmall, rlarge = [], [], [], [], []
    Yboth_, Ysingle_ = [], []
    delta_t = []
    taxa = list(data.df_microbiota_bool.columns)

    for i in range(m):
        if len(Yboth[i]) * len(Ysingle[i]) > 0:
            Usmall, psmall = mannwhitneyu(Yboth[i], Ysingle[i], alternative='less')
            Ularge, plarge = mannwhitneyu(Yboth[i], Ysingle[i], alternative='greater')

            delta_t.append(np.median(Yboth[i]) - np.median(Ysingle[i]))
            rsmall.append(1 - 2 * Usmall / (len(Yboth[i]) * len(Ysingle[i])))
            rlarge.append(1 - 2 * Ularge / (len(Yboth[i]) * len(Ysingle[i])))

            index.append(taxa[i])
            Psmall.append(psmall)
            Plarge.append(plarge)
            Yboth_.append(Yboth[i])
            Ysingle_.append(Ysingle[i])

    Psmall = np.array(Psmall)
    Plarge = np.array(Plarge)

    Psmall, Plarge = dict(zip(index, Psmall)), dict(zip(index, Plarge))
    Ysingle, Yboth = dict(zip(index, Ysingle_)), dict(zip(index, Yboth_))
    delta_t = dict(zip(index, delta_t))
    rsmall, rlarge = dict(zip(index, rsmall)), dict(zip(index, rlarge))

    return Psmall, Plarge, Ysingle, Yboth, rsmall, rlarge, delta_t


def FindSignificantTaxa(P, R, Delta, eps):
    '''Returns a DataFrame of taxa that survive Bonferroni correction at level eps.

    Inputs:
        * P (dict): raw p-values keyed by full taxonomy string
        * R (dict): rank-biserial correlation coefficients keyed by taxonomy string
        * Delta (dict): median contact-duration differences keyed by taxonomy string
        * eps (float): significance threshold after Bonferroni correction

    Outputs:
        * significant (DataFrame): columns taxon, taxon_full, p_value,
          correlation, delta_t; sorted by corrected p-value
    '''
    n_taxa = len(P.keys())
    print(f'The test was run on {n_taxa} taxa')

    taxa = np.array(list(P.keys()))
    p = np.array([P[x] * n_taxa for x in taxa if P[x] * n_taxa < eps])
    r = np.array([R[x] for x in taxa if P[x] * n_taxa < eps])
    delta = np.array([Delta[x] for x in taxa if P[x] * n_taxa < eps])

    taxa_full = np.array([x for x in P if P[x] * n_taxa < eps])
    taxa = np.array([getName(x) for x in taxa_full])

    idx = np.argsort(p)
    p, taxa, taxa_full, r, delta = p[idx], taxa[idx], taxa_full[idx], r[idx], delta[idx]
    significant = pd.DataFrame(columns=['taxon', 'taxon_full', 'p_value', 'correlation', 'delta_t'])
    significant.taxon = taxa
    significant.taxon_full = taxa_full
    significant.p_value = p
    significant.correlation = r
    significant.delta_t = delta

    return significant


def getName(name):
    '''Extracts the most specific informative taxonomic label from a
    semicolon-delimited taxonomy string.

    Uninformative labels (e.g. uncultured_bacterium, metagenome) are stripped
    before traversing from the finest to the coarsest taxonomic level.

    Inputs:
        * name (str): semicolon-delimited taxonomy string

    Returns:
        * str: the deepest non-ambiguous taxonomic label, or 'unknown'
    '''
    name = name.replace('uncultured_bacterium', '')
    name = name.replace('metagenome', '')
    name = name.replace('Unassigned', '')
    name = name.replace('RF39', '')
    name = name.replace('uncultured_organism', '')

    name_ = name.split(';')
    i = len(name_) - 1
    flag = 0

    while flag == 0:
        if i < 0:
            name = 'unknown'
            flag = 1
        else:
            if len(name_[i].split('__')) == 1:
                i = i - 1
            else:
                if len(name_[i].split('.')) > 1:
                    i = i - 1
                else:
                    flag = 1
                    name = name_[i]
    return name


def GetSim(df_net, g, measure):
    '''Bins contact weights logarithmically and returns the mean weight and
    mean microbiota similarity per bin.

    Use: weight_bins, sim_bins = GetSim(df_net, g, measure)

    Inputs:
        * df_net (pandas DataFrame): contact network with 'weight', 'jaccard',
          and 'TaxonomicJaccard' columns
        * g (int): bin frequency; controls the number of bins (between 1 and 100;
          larger values yield fewer bins)
        * measure (str): microbiota similarity column to aggregate;
          one of 'jaccard' or 'TaxonomicJaccard'

    Outputs:
        * weight_bins (list): mean log(weight + 1) per bin
        * sim_bins (list): mean microbiota similarity per bin
    '''
    perc = [np.percentile(np.log(df_net.weight + 1), x) for x in range(1, 100) if x % g == 0]
    df_net['log_w'] = np.sum([np.log(df_net.weight) >= p for p in perc], axis=0)
    all_weights = set(df_net.log_w)
    x, y = [], []

    for w in all_weights:
        idx = df_net.log_w == w
        x.append(np.log(df_net[idx].weight + 1).mean())
        y.append(df_net[idx][measure].mean())

    return x, y


def LinearFit(x, y):
    '''Returns the slope and intercept of an ordinary least-squares linear fit.

    Inputs:
        * x, y (array-like): input variables

    Outputs:
        * alpha (float): slope
        * c (float): intercept
    '''
    x, y = np.array(x), np.array(y)
    alpha = (np.mean(x * y) - np.mean(x) * np.mean(y)) / np.var(x)
    c = np.mean(y) - alpha * np.mean(x)
    return alpha, c


def ComputeROC(sim_train, sim_test, y_train, y_test, classifier):
    '''Trains a classifier and returns TPR and FPR on the test set.

    Inputs:
        * sim_train, sim_test (array): feature matrices
        * y_train, y_test (array): binary labels (1 = contact above threshold)
        * classifier (str): 'svm' (RBF kernel) or 'random_forest'

    Outputs:
        * TPR (float): true positive rate
        * FPR (float): false positive rate
    '''
    if classifier == 'svm':
        res = svm.SVC(kernel='rbf').fit(sim_train, y_train)
    if classifier == 'random_forest':
        res = RandomForestClassifier(max_depth=2, random_state=0).fit(sim_train, y_train)

    pred = res.predict(sim_test)

    P, N = np.sum(y_test == 1), np.sum(y_test == 0)
    if P == 0:
        TPR, FPR = 0, 0
    elif N == 0:
        TPR, FPR = 1, 1
    else:
        TPR = np.sum((y_test == 1) & (pred == 1)) / P
        FPR = np.sum((y_test == 0) & (pred == 1)) / N

    return TPR, FPR


def GetUmapEmbedding(D, n, nn):
    '''Computes a 2D UMAP embedding from a precomputed distance matrix and
    evaluates its quality via Spearman correlation with the original distances.

    Inputs:
        * D (array): n×n precomputed distance matrix
        * n (int): number of samples
        * nn (int): number of nearest neighbours for UMAP

    Outputs:
        * X (array): n×2 UMAP embedding coordinates
        * r (float): Spearman correlation between embedded and original distances
        * p (float): p-value of the Spearman correlation
    '''
    X = umap.UMAP(metric='precomputed', n_neighbors=nn).fit_transform(D)

    Dumap = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            Dumap[i, j] = np.linalg.norm(X[i] - X[j])

    res = spearmanr(Dumap.flatten(), D.flatten())
    return X, res.statistic, res.pvalue


def ComputeR2Permanova(res, D, data):
    '''Computes the R² effect size for a PERMANOVA result.

    R² is derived from the pseudo-F statistic as the fraction of total
    variance explained by the grouping factor.

    Inputs:
        * res: PERMANOVA result object (from skbio.stats.distance.permanova)
        * D (array): distance matrix used in the test
        * data (DATA): data container with df_microbiota (used to count groups)

    Outputs:
        * float: R² in [0, 1]
    '''
    F = res['test statistic']
    n, _ = D.shape
    g = len(data.df_microbiota['Person ID'].unique())
    return (1 + (n - g) / ((g - 1) * F)) ** (-1)


def getNameByHand(name, level):
    '''Extracts the taxonomic label at a specified rank from a
    semicolon-delimited taxonomy string.

    Inputs:
        * name (str): semicolon-delimited taxonomy string
        * level (str): taxonomic rank prefix; one of 'd' (domain), 'p' (phylum),
          'c' (class), 'g' (genus), 's' (species)

    Returns:
        * str: the label at the requested rank, stripped of trailing suffixes
    '''
    names = name.split(';')
    flag = 0
    i = len(names) - 1

    while flag == 0:
        x = names[i]
        x_split = x.split('__')
        if x_split[0] == level:
            flag = 1
            name_ = x
        else:
            i = i - 1

    return name_.split('.')[0]
