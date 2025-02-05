import numpy as np
from scipy.stats import mannwhitneyu
from collections import Counter
import matplotlib.pyplot as plt
from copy import copy
from scipy.special import betainc
from itertools import combinations
import pandas as pd
from sklearn.metrics import log_loss
from scipy.stats import spearmanr
from sklearn import svm


def distance(x, y, dist_type):

    if dist_type == 'jaccard':
        return jaccard_similarity(x,y)

    elif dist_type == 'js':
        return jensen_shannon(x,y)
    
    elif dist_type == 'bc':
        return bray_curtis(x,y)
    
    elif dist_type == 'cos':
        return cosine_similarity(x,y)
    
    elif dist_type == 'logloss':
        x_, y_ = np.sign(x), np.sign(y)
        return log_loss(x_, y_)
    
    elif dist_type == 'spearman':
        r, _ = spearmanr(x,y)
        return r


def jaccard_similarity(x1, x2):
    '''returns the jaccard similarity between two sets'''

    xx1 = (x1 != 0).astype(int)
    xx2 = (x2 != 0).astype(int)
    return (xx1 & xx2).sum() / (xx1 | xx2).sum()


def log0(x):

    y = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] > 0:
            y[i] = np.log(x[i])
        else:
            y[i] = 0

    return y


def dkl(x, y):
    '''computes the Kullback-Leibler divergence'''
    return np.sum(x*log0(x/y))

def cosine_similarity(x,y):
    '''computes the cosine similarity between x and y'''
    return x@y/np.sqrt((x@x)*(y@y))

def jensen_shannon(xx1, xx2):
    '''computes the jensen_shannon divergence'''
   
    if xx1.sum() > 0:
        xx1 = xx1 / xx1.sum()
    
    if xx2.sum() > 0:
        xx2 = xx2 / xx2.sum()

    xxm = (xx1 + xx2) / 2
    return 0.5 * dkl(xx1, xxm).sum() + 0.5 * dkl(xx2, xxm).sum()


def bray_curtis(x1, x2):
    '''computes the bray-curtis dissimilarity'''
    return 1.0 - 2.0 * np.min([x1,x2], axis=0).sum() / (x1.sum() + x2.sum())



def FindSignificantTaxa(data):
    '''This function finds the taxa for which we have a significant difference between the contact duration distribution of pairs both having that taxa and all others
    
    Use: Psmall, Plarge, Ysingle, Yboth = FindSignificantTaxa(data)
    
    Outputs:
        * Psmall: pvalues (one per bacterium) of rejection of the alternative 'less'
        * Plarge: pvalues (one per bacterium) of rejection of the alternative 'greater'
        * Yboth: list of arrays of the log contact duration of nodes both having that bacterium
        * Ysingle: list of arrays of the log contact duration of nodes pairs in which either of the two has the bacterium
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

    df_net = pd.DataFrame(np.array([idx1, idx2, weights]).T, columns = ['pid', 'pid2', 'weight'])
    df_net.weight = df_net.weight.astype(float)

    A = data.df_microbiota_bool.loc[df_net.pid].values
    B = data.df_microbiota_bool.loc[df_net.pid2].values
    IDX = (A * B) > 0
    _, m = IDX.shape

    Yboth = [df_net.weight.values[IDX[:,i]] for i in range(m)]
    Ysingle = [df_net.weight.values[~IDX[:,i]] for i in range(m)]



    Psmall, Plarge = [], []
    for i in range(m):
        if len(Yboth[i]) > 30 and len(Ysingle[i]) > 30:
            _, psmall = mannwhitneyu(Yboth[i], Ysingle[i], alternative = 'less')
            _, plarge = mannwhitneyu(Yboth[i], Ysingle[i], alternative = 'greater')

        else:
            psmall, plarge = 1, 1

        Psmall.append(psmall)
        Plarge.append(plarge)

    Psmall = np.array(Psmall)
    Plarge = np.array(Plarge)

    return Psmall, Plarge, Ysingle, Yboth
        

def ComputeROC(sim_train, sim_test, y_train, y_test):
    '''This function performs an SVM on (sim_train, y_train) and computes the TRP, FPR on the test set'''
    res = svm.SVC().fit(sim_train, y_train)
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