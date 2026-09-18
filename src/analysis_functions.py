from skbio.diversity import beta_diversity
import numpy as np
import pandas as pd
from itertools import combinations, permutations
from copy import copy
from scipy.stats import mannwhitneyu
from skbio.stats.distance import DistanceMatrix, permanova
import re
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm

class Model():
    def __init__(self, df_net, G, microbiota_all, microbiota, asv_ids, tree, days):
        self.df_net = df_net
        self.G = G
        self.microbiota_all = microbiota_all
        self.microbiota = microbiota
        self.asv_ids = asv_ids
        self.tree = tree
        self.days = days

def PrepareLevels(data, persistences = [1, 3], ths = [0, 1, 10]):

    Models = dict()

    for persistence in persistences:
        for th in ths:
            M = Model(None, None, None, None, None, None, None)
            M.df_net = data.df_net
            M.G = data.G
            M.tree = data.tree
            M.days = data.days

            microbiota_filtered, asv_ids_filtered = Filter(data, persistence, th)
            M.microbiota_all = microbiota_filtered
            M.asv_ids = asv_ids_filtered

            all_days = np.sort(np.unique([x[1] for x in data.microbiota_all.index]))
            M.microbiota = dict(zip(all_days, [microbiota_filtered.xs(day, level='Day code') for day in all_days]))

            Models[(persistence, th)] = M

            for metric in ['jaccard', 'braycurtis', 'weighted_unifrac', 'unweighted_unifrac']:
                Update_dfnet(Models[(persistence, th)], metric, column = metric)

    return Models


def Filter(data, persistence, th):

    microbiota_filtered = copy(data.microbiota_all)
    microbiota_filtered = FilterPersistentTaxa(microbiota_filtered, persistence)
    microbiota_filtered, asv_ids_filtered = FilterPrevalence(microbiota_filtered, data.asv_ids, th)

    return microbiota_filtered, asv_ids_filtered


def FilterPrevalence(microbiota_all, asv_ids, th):
    '''Drops ASVs (columns) present in at most `th` samples across the whole dataset. Singletons
    like this are a typical signature of sequencing/PCR noise.

    Use: microbiota_all, asv_ids = FilterPrevalence(microbiota_all, asv_ids, th)

    Inputs:
        * th (int): taxa appearing in `th` samples or fewer, across the whole dataset, are dropped.
    '''

    prevalence = (microbiota_all.values > 0).sum(axis=0)
    keep = prevalence > th

    microbiota_all = microbiota_all.loc[:, keep]
    asv_ids = np.asarray(asv_ids)[keep]

    return microbiota_all, asv_ids

def FilterPersistentTaxa(microbiota_all, persistence):
    '''Zeroes out, for each person individually, any taxon that recurs in fewer than
    `persistence` of that person's own first-wave samples (a transient detection is more likely
    noise than a stable member of that person's microbiota). Unlike FilterPrevalence, this is
    evaluated per person: a taxon can be kept for one person and dropped for another.

    Use: microbiota_all = FilterPersistentTaxa(microbiota_all, persistence)
    '''

    return microbiota_all * (np.sign(microbiota_all).groupby('Person ID').sum() >= persistence)



def Update_dfnet(data, dist_type, column=None):
    '''Rebuilds data.df_net as a complete pid/pid2 table (adding zero-weight edges for pairs that
    never interacted) and adds a column with, for each pair, their microbiota distance averaged
    across their first-wave (up to 5) measurements. Can be called repeatedly with different
    dist_type values to add several distance columns.

    Use: Update_dfnet(data, dist_type, column=None)

    Inputs:
        * data (LoadData.ReturnValue): loaded data, must have df_net, G, microbiota_all, asv_ids, tree
        * dist_type (str): skbio beta diversity metric, e.g. 'jaccard', 'braycurtis',
          'weighted_unifrac', or 'unweighted_unifrac'
        * column (str, optional): name of the column to store the result in data.df_net; defaults to dist_type
    '''

    column = column or dist_type

    # create a complete network with w = 0 for the missing edges
    idx1, idx2, weight = [], [], []
    all_nodes = np.unique(data.df_net[['pid', 'pid2']])
    for pid, pid2 in combinations(all_nodes, 2):
        idx1.append(pid)
        idx2.append(pid2)
        if tuple([pid, pid2]) in data.G.edges:
            weight.append(data.G[pid][pid2]['weight'])
        else:
            weight.append(0)

    df_net = pd.DataFrame({'pid': idx1, 'pid2': idx2, 'weight': weight})

    ids = data.microbiota_all.index.map(lambda t: '_'.join(map(str, t)))
    counts = data.microbiota_all.values

    kwargs = {}
    if dist_type in ('weighted_unifrac', 'unweighted_unifrac'):
        kwargs = dict(taxa=data.asv_ids, tree=data.tree)

    elif dist_type == 'braycurtis':
        counts = counts / counts.sum(axis=1, keepdims=True)

    dm = beta_diversity(metric=dist_type, counts=counts, ids=ids, validate=True, **kwargs)

    df_net[column] = df_net.apply(lambda x: average_distance(dm, x.pid, x.pid2, data.days), axis = 1)

    # preserve similarity columns computed by previous calls with other dist_type values
    # (excluding `column` itself, so re-running with the same dist_type overwrites instead of duplicating)
    prev_cols = [c for c in data.df_net.columns if c not in ('pid', 'pid2', 'weight', column)]
    if prev_cols:
        df_net = df_net.merge(data.df_net[['pid', 'pid2'] + prev_cols], on = ['pid', 'pid2'], how = 'left')

    data.df_net = df_net

    return

def average_distance(dm, pid, pid2, days):
        dists = []
        for day in days:
            a, b = f'{pid}_{day}', f'{pid2}_{day}'
            if a in dm.ids and b in dm.ids:
                dists.append(dm[a, b])
        return np.mean(dists)

def People_vs_time(data, dist_type):
    '''Computes, for every pair of samples, the microbiota distance between repeated measurements
    of the same person (sim) versus between different people (nsim) - the basis for testing
    whether an individual's microbiota is more stable across time than across people.

    Use: sim, nsim = People_vs_time(data, dist_type)

    Inputs:
        * data (LoadData.ReturnValue): loaded data, must have microbiota_all, asv_ids, tree
        * dist_type (str): skbio beta diversity metric, e.g. 'jaccard', 'braycurtis',
          'weighted_unifrac'.

    Outputs:
        * sim (list): pairwise distances between repeated samples of the same person
        * nsim (list): pairwise distances between samples of different people
    '''

    sim, nsim = [], []
    ids = data.microbiota_all.index.map(lambda t: '_'.join(map(str, t)))
    counts = data.microbiota_all.values

    kwargs = {}
    if dist_type in ('weighted_unifrac', 'unweighted_unifrac'):
        kwargs = dict(taxa=data.asv_ids, tree=data.tree)
    elif dist_type == 'braycurtis':
        # bray-curtis is sensitive to sequencing depth, so normalize to relative abundances first
        counts = counts / counts.sum(axis=1, keepdims=True)

    dm = beta_diversity(metric=dist_type, counts=counts, ids=ids, validate=True, **kwargs)

    for a in ids:
        for b in ids:
            if a != b:
                if a.split('_')[0] == b.split('_')[0]:
                    sim.append(dm[a, b])
                else:
                    nsim.append(dm[a, b])

    return sim, nsim

def df_net_strong_ties(data, column, n_sim=1500):
    '''Null model for the mean microbiota similarity among pairs whose contact duration exceeds
    increasing thresholds: contact durations are randomly reshuffled across pairs `n_sim` times,
    breaking any true link between contact duration and similarity, so the resulting distribution
    shows what the threshold curve would look like by chance alone.

    Use: Sim, tv = df_net_strong_ties(data, column)

    Inputs:
        * data (LoadData.ReturnValue or Model): must have df_net with a weight column and the
          requested similarity column
        * column (str): name of the data.df_net column to use, e.g. 'jaccard', 'braycurtis'
        * n_sim (int): number of randomizations

    Outputs:
        * Sim (list): for each randomization, the mean of column at each threshold in tv
        * tv (array): threshold values (in hours)
    '''

    # set of thresholds considered (the unit is hours)
    tv = np.linspace(0, 6, 6)

    # randomize and store in X the mean similarity
    df_net_rdn = copy(data.df_net)
    x = copy(df_net_rdn.weight.values)
    Sim = []

    for i in range(n_sim):
        print(i, end = '\r')

        np.random.shuffle(x)
        df_net_rdn.weight = x
        Sim.append([df_net_rdn[df_net_rdn.weight*10/3600 > t][column].mean() for t in tv])

    return Sim, tv


def GetMicrobiotaBool(data):
    '''Builds a (person x taxon) boolean presence table: 1 if a person ever presented a taxon
    across their first-wave (up to 5) measurements, 0 otherwise

    Use: df_microbiota_bool = GetMicrobiotaBool(data)

    Inputs:
        * data (LoadData.ReturnValue): loaded data, must have microbiota_all, G

    Outputs:
        * df_microbiota_bool (DataFrame): index = Person ID (as str, to match data.G's node dtype),
          columns = taxa, values in {0, 1}
    '''

    return np.sign(data.microbiota_all).groupby('Person ID').sum()


def GetTaxaPValues(data):
    '''For each taxon, tests whether pairs who both carry it have significantly different contact
    durations than pairs where only one of them does (Mann-Whitney U, both tails).

    Use: Psmall, Plarge, Ysingle, Yboth, rsmall, rlarge, delta_t = GetTaxaPValues(data)

    Inputs:
        * data (LoadData.ReturnValue): loaded data, must have df_net, G, and df_microbiota_bool
          (set data.df_microbiota_bool = GetMicrobiotaBool(data) beforehand)

    Outputs:
        * Psmall, Plarge (dict): p-values for the 'less'/'greater' alternatives, per taxon
        * rsmall, rlarge (dict): rank-biserial effect size matching Psmall/Plarge, per taxon
        * Yboth (dict): contact durations of pairs where both members carry the taxon
        * Ysingle (dict): contact durations of pairs where exactly one member carries the taxon
        * delta_t (dict): median(Yboth) - median(Ysingle), per taxon
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

    df_net = pd.DataFrame({'pid': idx1, 'pid2': idx2, 'weight': weights})

    A = data.df_microbiota_bool.loc[df_net.pid].values
    B = data.df_microbiota_bool.loc[df_net.pid2].values
    IDXboth = (A * B) > 0
    IDXone = (np.sign(A) + np.sign(B)) * (1 - np.sign(A) * np.sign(B)) > 0
    _, m = IDXboth.shape

    Yboth = [df_net.weight.values[IDXboth[:,i]] for i in range(m)]
    Ysingle = [df_net.weight.values[IDXone[:,i]] for i in range(m)]

    index, Psmall, Plarge, rsmall, rlarge = [], [], [], [], []
    Yboth_, Ysingle_ = [], []
    delta_t = []
    taxa = list(data.df_microbiota_bool.columns)

    for i in range(m):   
        if len(Yboth[i]) * len(Ysingle[i]) > 0:
            Usmall, psmall = mannwhitneyu(Yboth[i], Ysingle[i], alternative = 'less')
            Ularge, plarge = mannwhitneyu(Yboth[i], Ysingle[i], alternative = 'greater')

            delta_t.append(np.median(Yboth[i]) - np.median(Ysingle[i]))
            rsmall.append(1 - 2*Usmall/(len(Yboth[i])*len(Ysingle[i])))
            rlarge.append(1 - 2*Ularge/(len(Yboth[i])*len(Ysingle[i])))

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


def FindSignificantTaxa(P, R, Delta, eps = 0.05):
    '''Returns a DataFrame of the taxa whose p-value in P survives Bonferroni correction at level eps'''

    n_taxa = len(P.keys())
    print(f'The test was run on {n_taxa} taxa')

    taxa = np.array(list(P.keys()))
    p = np.array([P[x]*n_taxa for x in taxa if P[x]*n_taxa < eps])
    r = np.array([R[x] for x in taxa if P[x]*n_taxa < eps])
    delta = np.array([Delta[x] for x in taxa if P[x]*n_taxa < eps])
    
    taxa_full = np.array([x for x in P if P[x]*n_taxa < eps])
    taxa =  np.array([getName(x) for x in taxa_full])
    
    idx = np.argsort(p)
    p, taxa, taxa_full, r, delta = p[idx], taxa[idx], taxa_full[idx], r[idx], delta[idx]
    significant = pd.DataFrame(columns = ['taxon', 'taxon_full', 'p_value', 'correlation', 'delta_t'])
    significant.taxon = taxa
    significant.taxon_full = taxa_full
    significant.p_value = p
    significant.correlation = r
    significant.delta_t = delta

    return significant


def getName(name):
    '''Returns the last available taxonomic level of a semicolon-separated lineage string
    (e.g. "d__Bacteria;...;g__Streptococcus;s__.11"). Each level is stripped only of its
    trailing per-ASV disambiguator (e.g. ".11", ".1") before being checked for content, so a
    genuine species-level epithet (e.g. "Corynebacterium_durum.10" -> "Corynebacterium_durum")
    is kept, and only a level that is truly blank or a known placeholder after that falls back
    to the next shallower level.
    '''

    name = name.replace('uncultured_bacterium', '')
    name = name.replace('metagenome', '')
    name = name.replace('Unassigned', '')
    name = name.replace('RF39', '')
    name = name.replace('uncultured_organism', '')

    name_ = name.split(';')
    i = len(name_)-1

    flag = 0

    while flag == 0:
        if i < 0:
            name = 'unknown'
            flag = 1
        else:
            level = name_[i]
            if len(level.split('__')) == 1:
                i = i-1
            else:
                rank, _, value = level.partition('__')
                value = re.sub(r'\.\d+$', '', value)
                if value == '':
                    i = i-1
                else:
                    flag = 1
                    name = f'{rank}__{value}'
    return name


def getNameByHand(name, level):
    '''Shortens a lineage string to its epithet at the given rank, e.g.
    getNameByHand('...;g__Streptococcus;...', 'g') -> 'g_Streptococcus' '''

    epithet = name.split(f'{level}__')[1].split(';')[0].split('_')[0]

    return f'{level}_{epithet}'


def LinearFit(x, y):
    '''Ordinary least-squares fit of y = alpha * x + c'''

    x, y = np.array(x), np.array(y)
    alpha = (np.mean(x*y) - np.mean(x)*np.mean(y)) / np.var(x)
    c = np.mean(y) - alpha*np.mean(x)

    return alpha, c


def ComputeROC(sim_train, sim_test, y_train, y_test, classifier):
    '''Fits `classifier` on (sim_train, y_train) and returns its TPR, FPR on the test set'''

    if classifier == 'svm':
        res = svm.SVC(kernel = 'rbf').fit(sim_train, y_train)    
    if classifier == 'random_forest':
        res = RandomForestClassifier(max_depth=2, random_state=0).fit(sim_train, y_train)
    
    pred = res.predict(sim_test)

    P, N = np.sum(y_test == 1), np.sum(y_test == 0)
    if P == 0:
        TPR, FPR = 0, 0
    elif N == 0:
        TPR, FPR = 1, 1
    else:
        TPR = np.sum((y_test == 1) & (pred == 1))/P
        FPR = np.sum((y_test == 0) & (pred == 1))/N

    return TPR, FPR


def BuildDistanceMatrix(data, dist):

    ids = data.microbiota_all.index.map(lambda t: '_'.join(map(str, t)))
    counts = data.microbiota_all.values

    kwargs = {}
    if dist in ('weighted_unifrac', 'unweighted_unifrac'):
        kwargs = dict(taxa=data.asv_ids, tree=data.tree)
    elif dist == 'braycurtis':
        # bray-curtis is sensitive to sequencing depth, so normalize to relative abundances first
        counts = counts / counts.sum(axis=1, keepdims=True)

    dm = beta_diversity(metric=dist, counts=counts, ids=ids, validate=True, **kwargs)

    return dm.data


def RunPermanova(data, dist):
    '''Runs PERMANOVA on the microbiota distance matrix twice: grouping samples by sampling day,
    and grouping samples by person. A low R2/p-value for "day" and a high one for "person" would
    indicate that microbiota composition is driven more by individual identity than by when the
    sample was taken.

    Use: R2_day, p_day, R2_person, p_person = RunPermanova(data, dist)
    '''

    ids = data.microbiota_all.index.map(lambda t: '_'.join(map(str, t)))
    D = BuildDistanceMatrix(data, dist)

    day_codes = data.microbiota_all.index.get_level_values('Day code')
    person_id_grouping = data.microbiota_all.index.get_level_values('Person ID')

    D_arr = np.asarray(D)
    perm_dm = DistanceMatrix(D_arr, ids=list(ids))
    res1 = permanova(perm_dm, grouping=day_codes)
    res2 = permanova(perm_dm, grouping=person_id_grouping)

    return ComputeR2Permanova(res1, D, day_codes), res1['p-value'], ComputeR2Permanova(res2, D, person_id_grouping), res2['p-value']


def ComputeR2Permanova(res, D, grouping):
    '''Converts a PERMANOVA pseudo-F statistic into an R2 effect size'''

    F = res['test statistic']
    n, _ = D.shape
    g = len(np.unique(grouping))

    return (1 + (n-g)/((g-1)*F))**(-1)
