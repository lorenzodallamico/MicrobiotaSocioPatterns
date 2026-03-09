import numpy as np


def jaccard_similarity(x1, x2):
    '''Returns the Jaccard similarity between two binary vectors.

    Inputs:
        * x1, x2 (array): abundance or binary presence/absence vectors

    Returns:
        * float: Jaccard similarity (intersection over union of non-zero elements)
    '''
    xx1 = (x1 != 0).astype(int)
    xx2 = (x2 != 0).astype(int)
    return (xx1 & xx2).sum() / (xx1 | xx2).sum()


def bray_curtis(x1, x2):
    '''Returns the Bray-Curtis dissimilarity between two abundance vectors.

    Inputs:
        * x1, x2 (array): non-negative abundance vectors

    Returns:
        * float: Bray-Curtis dissimilarity in [0, 1]
    '''
    return 1.0 - 2.0 * np.min([x1, x2], axis=0).sum() / (x1.sum() + x2.sum())


def TaxonomicJaccard(D, i, j):
    '''Returns the unweighted Taxonomic Jaccard distance between two microbiota samples.

    The distance is computed as one minus the Jaccard similarity of the sets
    of all taxonomic labels (across all levels) present in each sample.

    Inputs:
        * D (dict): branch-set dictionary produced by getBranches()
        * i, j: sample indices (keys in D)

    Returns:
        * float: Taxonomic Jaccard distance in [0, 1]
    '''
    X, Y = D[i], D[j]
    return 1 - len(X.intersection(Y)) / len(X.union(Y))


def getBranches(df):
    '''Builds the branch-set dictionary required for computing Taxonomic Jaccard distances.

    For each sample, collects all individual taxonomic labels present across
    every level of the semicolon-delimited taxonomy strings.

    Inputs:
        * df (pandas DataFrame): abundance table; index = sample IDs,
          columns = semicolon-delimited taxonomy strings

    Returns:
        * D (dict): maps each sample index to the set of all taxonomic labels
          present in that sample
    '''
    D = dict()

    for i in df.index:
        x = df.loc[i]
        x = x[x > 0]
        taxa = list(x[x > 0].index)
        X = []

        for tp in taxa:
            t = tp.split(';')
            for t_ in t:
                X.append(t_)

        D[i] = set(X)

    return D
