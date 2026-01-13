import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))
import traversal as trv
import bfutils as bfu
import graphkit as gk
import zickle as zkl
import numpy as np
import itertools as iter
import statsmodels.api as sm
import linear_model as lm
import ipdb
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore, norm

def independent(y,X,pval=0.05):
    X  = sm.add_constant(X)
    est  = sm.OLS(y,X).fit()
    return est.pvalues[1] > pval

def kernel(z):
    if np.abs(z) > 1.: return 0.
    return .5

def residuals_(x,y,z):
    N = len(x)
    pd = map(kernel,pdist(z.T))
    PD = squareform(pd)
    sumsx = np.dot(PD,x) + 0.5*x
    sumsy = np.dot(PD,y) + 0.5*y
    weights = np.sum(PD, axis=1) + 0.5
    residualsx = x - sumsx/weights
    residualsy = y - sumsy/weights

    return residualsx, residualsy

def moment22(x,y): return np.dot(x*x,y*y)/len(x)

def fdr(alpha, pvalues):
    m = len(pvalues)
    c = np.cumsum(1./(np.arange(m)+1))
    pcompare = (pvalues <= alpha*(np.arange(1,m+1)/(c*(m+1))))
    idx = np.where(pcompare == True)[0]
    if len(idx) == 0: return -1
    return idx[-1]

def fdrQ(alpha,pvalues):
    pvalues = np.sort(pvalues)
    pvalues = pvalues[~np.isnan(pvalues)]
    min  = np.nan if len(pvalues) == 0 else pvalues[0]
    high = 1.0
    low  = 0.
    q    = alpha
    while (high - low) > 1e-5:
        midpoint = (high + low)/2.0
        q = midpoint
        cutoff = pvalues[fdr(q, pvalues)]

        if cutoff < min:
            low = midpoint
        elif cutoff > min:
            high = midpoint
        else:
            low  = midpoint
            high = midpoint
    return q

def fdrCutoff(alpha, pvalues):
    pvalues = np.sort(pvalues)
    k = fdr(alpha,pvalues)
    if k < 0: return 0
    return pvalues[k]

def np_fisherZ(x,y,r):
    z = 0.5 * (np.log(1.0 + r) - np.log(1.0 - r))
    w = np.sqrt(len(x)) * z
    x_ = zscore(x)
    y_ = zscore(y)
    t2 = moment22(x_,y_)
    t = np.sqrt(t2)
    p = 2. * (1. - norm.cdf(np.abs(w), 0.0, t))
    return p

def independent_(x,y, alpha = 0.05):
    x = x[~np.isnan(x)]
    y = y[~np.isnan(y)]

    # For PC, should not remove the edge for this reason.
    if len(x) < 10: return False

    ps = []
    for i in range(15):
        for j in range(15):
            x_ = x**i
            y_ = y**j
            r = np.corrcoef(x_,y_)[0,1]
            #r = max(min(r,1),-1) # Tetrad had this
            p = np_fisherZ(x_,y_,r)
            #if not np.isnan(p):
            ps.append(p)

    if not ps: return True
    return fdrCutoff(alpha, ps) > alpha


def addallb(g):
    n = len(g)
    for i in range(n):
        for j in range(n):
            if str(j+1) in g[str(i+1)]:
                g[str(i+1)][str(j+1)].add((2,0))
            else:
                g[str(i+1)][str(j+1)] = set([(2,0)])
    return g

def dpc(data, pval=0.05):
    n = data.shape[0]
    # stack the data: first n rows is t-1 slice, the next n are slice t
    data = np.asarray(np.r_[data[:,:-1],data[:,1:]])

    def tetrad_cind_(y,x,condset=[], alpha=0.01, shift=0):
        y = data[n+int(y)-1,:]
        x = data[shift+int(x)-1,:]
        if condset:
            X  = data[condset,:]
            ry, rx = residuals_(y,x,X)
        else:
            ry, rx = [y,x]
        return independent_(ry, rx, alpha = alpha)

    def cind_(y,x, condset=[], pval=pval, shift=0):
        yd = data[n+int(y)-1,:].T
        X  = data[[shift+int(x)-1]+condset,:].T
        return independent(yd, X, pval=pval)

    def cindependent(y, x, counter, parents=[], pval=pval):
        for S in [j for j in iter.combinations(parents,counter)]:
            if cind_(y, x, condset=list(S), pval=pval): return True
            #if tetrad_cind_(x, y, condset=list(S), alpha=pval): return True                
        return False

    def bindependent(y, x, parents=[], pval=pval):
        return cind_(y, x, condset=parents, pval=pval, shift=n)
        #return tetrad_cind_(y, x, condset=parents, alpha=pval, shift=n)

    def prune(elist, mask, g):
        for e in mask:
            g[e[0]][e[1]].remove((0,1))
            elist.remove(e)
        gk.clean_leaf_nodes(g)

    g  = gk.superclique(n)
    gtr= bfu.gtranspose(g)

    el = gk.edgelist(g)
    for counter in range(n):
        to_remove = []
        for e in el:
            ppp = [int(k)-1 for k in gtr[e[1]] if k != e[0]]
            if counter <= len(ppp):
                if cindependent(e[1], e[0], counter, parents=ppp, pval=pval):
                    to_remove.append(e)
                    gtr[e[1]].pop(e[0],None)
        prune(el, to_remove, g)

    bel = [map(lambda k: str(k+1), x) for x in iter.combinations(range(n),2)]
    for e in bel:
        ppp = list(set(gtr[e[0]].keys()) | set(gtr[e[1]].keys()))
        ppp = map(lambda x: int(x)-1, ppp)
        if bindependent(e[0], e[1], parents=ppp, pval=pval):
            g[e[0]][e[1]].remove((2,0))
            g[e[1]][e[0]].remove((2,0))
    gk.clean_leaf_nodes(g)

    return g


# Local Variables:
# mode: python
# python-indent-offset: 4
# End:
