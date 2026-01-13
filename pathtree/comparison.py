import ecj
import scipy
import numpy
import operator
import networkx as nx
#from progressbar import ProgressBar, Percentage
numpy.random.RandomState()
import bfutils as bfu
import numpy as np
import gmpy2 as gmp

def num2CG(num,n):
    """num2CG - converts a number  whose binary representaion encodes edge
    presence/absence into a compressed graph representaion

    """
    n2 = n*n
    G = {'%i'%(i+1):{} for i in range(n)}
    if num == 0: return G
    bl = gmp.bit_length(num)
    idx = [n2-i-1 for i in range(bl) if num & (1<<i)]
    idx = np.unravel_index(idx,(n,n))
    x = idx[0]+1
    y = idx[1]+1
    for i in range(len(x)):
        G['%i' % x[i]]['%i' % y[i]] = set([(0,1)])
    return G

def hasSelfLoops(G):
    for u in G:
        if G[u].has_key(u):
            return True
    return False

def randSCC(n):
    G = num2CG(scipy.random.randint(2**(n**2)),n)
    while (len(ecj.scc(G)) > 1) or gcd4scc(G)>1:
        G = num2CG(scipy.random.randint(2**(n**2)),n)
    return G

def SM_fixed(Gstar,G, iter=5):
    compat = []
    for j in range(0,iter):
        if Gstar == ecj.undersample(G,j):
            compat.append(j)
    return compat
def SM_converging(Gstar,G):
    """Gstar is the undersampled reference graph, while G is the starting
    graph. The  code searches  over all undersampled  version of  G to
    find all matches with Gstar
    """
    compat = []
    GG = G
    Gprev = G
    if G == Gstar: return [0]
    j = 1
    G = ecj.undersample(GG,j)
    while not (G == Gprev):
        if Gstar == G: compat.append(j)
        j += 1
        Gprev = G
        G = ecj.undersample(GG,j)
    return compat
def searchMatch(Gstar,G, iter=5):
    if gcd4scc(G) >1: return SM_fixed(Gstar, G, iter=iter)
    return SM_converging(Gstar, G)
def hasSink(G):
    return not reduce(operator.and_, [bool(G[n]) for n in G], True)
def hasRoot(G): return hasSink(ecj.tr(G))
def isSclique(G):
    n = len(G)
    for v in G:
        if sum([(0,1) in G[v][w] for w in G[v]]) < n: return False
        if sum([(2,0) in G[v][w] for w in G[v]]) < n-1: return False
    return True

def graph2nx(G):
    g = nx.DiGraph()
    for v in G:
        g.add_edges_from([(v,x) for x in G[v] if (0,1) in G[v][x]])
    return g
def nx2graph(G):
    g = {str(n+1):{} for n in G}
    for n in G:
        g['%i' % (n+1)] = {'%i' % (x+1):set([(0,1)]) for x in G[n]}
    return g

def gcd4scc(SCC):
    g = graph2nx(SCC)
    return ecj.listgcd(map(lambda x: len(x)-1, nx.simple_cycles(g)))

def compatibleAtU(uGstar):
    compat = []
    n = len(uGstar)
    numG = 2**(n**2)
    #pbar = Percentage()
    for i in range(1,numG):
        G = num2CG(i,n)
        #pbar.update(i+1)
        if len(ecj.scc(G)) > 1: continue
        l = searchMatch(uGstar,G, iter = 5)
        if l: compat.append((l,G))
    #pbar.finish()
    return compat
