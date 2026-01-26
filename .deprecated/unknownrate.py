#BFS implementation of subgraph and supergraph
#Gu to G1 algorithm
from itertools import combinations, permutations
from functools import wraps
import copy
import time
import sys,os
import numpy as np
import ipdb
import operator
from scipy.misc import comb
import math
import gmpy as gmp
import gmpy as gmp
from scipy.misc import comb
import zickle as zkl
import simpleloops as sls
import math
import load_loops
from matplotlib.cbook import flatten
from progressbar import ProgressBar, Percentage, \
    Bar, RotatingMarker, ETA, FileTransferSpeed

TOOLSPATH='./tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))


circp = zkl.load('circular_p.zkl')
alloops = load_loops.alloops

import pprint
import bfutils as bfu
import traversal as trv
import graphkit as gk
import comparison as cmp
import simpleloops as sl

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = trv.gsig(args[0])         # Signature: just the g
        #s = tool.signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def prune_conflicts(H, g, elist):
    """checks if adding an edge from the list to graph g causes a
    conflict with respect to H and if it does removes the edge
    from the list

    Arguments:
    - `H`: the undersampled graph
    - `g`: a graph under construction
    - `elist`: list of edges to check
    """
    l  = []
    for e in elist:
        gk.addanedge(g,e)
        if not bfu.call_u_conflicts(g, H): l.append(e)
        gk.delanedge(g,e)
    return l

def eqclass(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    @memo
    def addedges(g,H,edges):
        if edges:
            nedges = prune_conflicts(H, g, edges)
            n = len(nedges)

            if n == 0: return None

            for i in range(n):
                gk.addanedge(g,nedges[i])
                if bfu.call_u_equals(g, H): s.add(bfu.g2num(g))
                addedges(g,H,nedges[:i]+nedges[i+1:])
                gk.delanedge(g,nedges[i])
    edges = gk.edgelist(gk.complement(g))
    addedges(g,H,edges)

    return s-set([None])

# these two functions come fromt his answer:
# http://stackoverflow.com/a/12174125
def set_bit(value, bit): return value | (1<<bit)
def clear_bit(value, bit): return value & ~(1<<bit)
def e2num(e,n): return (1<<(n*n +n - int(e[0],10)*n - int(e[1],10)))
def le2num(elist,n):
    num = 0
    for e in elist:
        num |= e2num(e,n)
    return num
def ekey2e(ekey,n):
    idx = np.unravel_index(n*n - ekey .bit_length() - 1 + 1,(n,n))
    idx = tuple([x+1 for x in idx])
    return ('%i %i'%idx).split(' ')

def cacheconflicts(num, cache):
    """Given a number representation of a graph and an iterable of
    conflicting subgraphs return True if the graph conflicts with any
    of them and false otherwise

    Arguments:
    - `num`: the number representation of a graph
    - `cache`: an iterable of number representations of conflicting
      graphs
    """
    conflict = False
    for c in cache:
        if num & c == c:
            return True
    return False

class nobar:
    def update(self,c): return None
    def finish(self): return None

def start_progress_bar(iter, n, verbose = True):
    if verbose:
        pbar = ProgressBar(widgets=['%3s' % str(iter) +
                                '%10s' % str(n)+' ',
                                Bar('-'), ' '],
                        maxval=n).start()
    else:
        pbar = nobar()
    return pbar

def add2set_loop(ds, H, cp, ccf, iter=1, verbose=True,
                 capsize=100, currsize=0):
    n = len(H)
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds), verbose = verbose)

    c = 0

    for gnum in ds:
        c += 1
        pbar.update(c)
        gset = set()
        eset = set()
        for sloop in ds[gnum]:
            if sloop & gnum == sloop: continue
            num = sloop | gnum
            if sloop in ccf and skip_conflictors(num,ccf[sloop]):continue
            if not num in s:
                g = bfu.num2CG(num, n)
                if not bfu.call_u_conflicts(g, H):
                    s.add(num)
                    gset.add((num,sloop))
                    eset.add(sloop)
                    if bfu.call_u_equals(g, H):
                        ss.add(num)
                        if capsize <= len(ss)+currsize: return dsr, ss

        #for gn,e in gset:
        #   if e in cp:
        #       dsr[gn] = eset - cp[e] - set([e])
        #   else:
        #       dsr[gn] = eset - set([e])
        for gn in gset: dsr[gn[0]] = eset - set([gn[1]])
    pbar.finish()
    return dsr, ss

def add2set_(ds, H, cp, ccf, iter=1, verbose=True, capsize=100):
    n = len(H)
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds), verbose = verbose)

    c = 0

    for gnum in ds:
        g = bfu.num2CG(gnum, n)
        c += 1
        pbar.update(c)
        glist = []
        elist = []
        eset = set()
        for e in ds[gnum]:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                ekey = (1<<(n2 - int(e[0],10)*n - int(e[1],10)))
                if ekey in ccf and skip_conflictors(num,ccf[ekey]):
                    gk.delanedge(g,e)
                    continue
                if not num in s:
                    s.add(num)
                    if not bfu.call_u_conflicts(g, H):
                        #cf, gl2 = bfu.call_u_conflicts2(g, H)
                        #if not cf:
                        glist.append((num,ekey))
                        elist.append(e)
                        eset.add(ekey)
                        if bfu.call_u_equals(g, H):
                            ss.add(num)
                            #if bfu.call_u_equals2(g, gl2, H): ss.add(num)
                        if capsize <= len(ss): break
                gk.delanedge(g,e)

        for gn,e in glist:
            if e in cp:
                dsr[gn] = [ekey2e(k,n) for k in eset - cp[e]]
            else:
                dsr[gn] = elist
        if capsize <= len(ss): return dsr, ss

    pbar.finish()
    return dsr, ss

def skip_conflictors(gnum, ccf):
    pss = False
    for xx in ccf:
        if xx&gnum == xx:
            pss = True
            break
    return pss

def bconflictor(e,H):
    n = len(H)
    s = set()
    for v in H:
        s.add(e2num((v,e[0]),n)|e2num((v,e[1]),n))
    return s

def conflictor(e,H):
    n = len(H)
    def pairs(e):
        ekey = e2num(e,n)
        return [ekey|e2num((e[0],e[0]),n),
                ekey|e2num((e[1],e[1]),n)]

    def trios(e,H):
        s = set()
        for v in H:
            if not v in e:
                s.add(e2num((e[0],v),n)|
                      e2num((v,e[1]),n)|
                      e2num((e[1],e[1]),n))
                s.add(e2num((e[0],e[0]),n)|
                      e2num((e[0],v),n)|
                      e2num((v,e[1]),n))
                s.add(e2num((e[0],v),n)|
                      e2num((v,v),n)|
                      e2num((v,e[1]),n))
        return s

    return trios(e,H).union(pairs(e))

def conflictor_set(H):
    s = set()
    for x in gk.inedgelist(H):  s = s | conflictor(x,H)
    for x in gk.inbedgelist(H): s = s | bconflictor(x,H)
    return s

def conflictors(H):
    s = conflictor_set(H)
    ds = {}
    num = reduce(operator.or_,s)
    for i in xrange(gmp.bit_length(num)):
        if num & 1<<i:
            ds[1<<i] = [x for x in s if x&(1<<i)]
    return ds

def may_be_true_selfloop(n,H):
    for v in H[n]:
        if v == n: continue
        if (0,1) in H[n][v] and not ((2,0) in H[n][v]): return False
    return True

def issingleloop(num):
    bl = gmp.bit_length(num)
    idx = [1 for i in xrange(bl) if num & (1<<i)]
    return len(idx) == 1

def nonbarren(H):
    for v in H:
        if H[v]: return v
    return False

def prune_loops(loops, H):
    l = []
    n = len(H)
    for loop in loops:
        g = bfu.num2CG(loop, n)
        x = [k for k in g if g[k]]
        if len(x) == 1:
            s = reduce(lambda x, s: s.union(x),
                       [H[x[0]][w] for w in H[x[0]]])
            if not (2,0) in s: continue
        if not bfu.call_u_conflicts_d(g, H): l.append(loop)        
    return l

def lconflictors(H, sloops=None):
    if not sloops: sloops = prune_loops(allsloops(len(H)),H)
    s = conflictor_set(H)
    ds = {}
    num = reduce(operator.or_,s)
    for i in xrange(gmp.bit_length(num)):
        if num & 1<<i:
            cset = [x for x in s if x&(1<<i)]
            for sloop in sloops:
                if sloop & 1<<i:
                    ds.setdefault(sloop,[]).extend(cset)
    return ds

def confpairs(H):
    n = len(H)
    g = {n:{} for n in H}
    d = {}

    edges = gk.edgelist(gk.complement(g))
    edges = prune_conflicts(H, g, edges)

    for p in combinations(edges,2):
        gk.addedges(g,p)
        if bfu.call_u_conflicts(g, H):
            n1 = e2num(p[0],n)
            n2 = e2num(p[1],n)
            d.setdefault(n1,set()).add(n2)
            d.setdefault(n2,set()).add(n1)
        gk.deledges(g,p)

    return d

def lconfpairs(H, cap=10, sloops=None):
    n = len(H)
    d = {}
    if not sloops: sloops = prune_loops(allsloops(len(H)),H)
    c = 0
    for p in combinations(sloops,2):
        g = bfu.num2CG(p[0]|p[1], n)
        if bfu.call_u_conflicts(g, H):
            d.setdefault(p[0],set()).add(p[1])
            d.setdefault(p[1],set()).add(p[0])
        if c >= cap: break
        c +=1
    return d


def iteqclass(H, verbose=True, capsize=100):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    if cmp.isSclique(H):
        print 'not running on superclique'
        return None
    g = {n:{} for n in H}
    s = set()
    Hnum = bfu.ug2num(H)
    if Hnum[1]==0: s.add(Hnum[0])

    cp = confpairs(H)
    ccf = conflictors(H)

    edges = gk.edgelist(gk.complement(g))
    ds = {bfu.g2num(g): edges}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    for i in range(len(H)**2):
        ds, ss = add2set_(ds, H, cp, ccf, iter=i,
                            verbose=verbose,
                            capsize=capsize)
        s = s | ss
        if capsize <= len(ss): break
        if not ds: break

    return s

def liteqclass(H, verbose=True, capsize=100, asl=None):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set([-1])
    g = {n:{} for n in H}
    s = set()



    if asl:
        sloops = asl
    else:
        sloops = prune_loops(allsloops(len(H)),H)

    cp  = []#lconfpairs(H, sloops=sloops)
    ccf = lconflictors(H, sloops=sloops)
    ds = {0: sloops}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    i=0
    while ds:
        ds, ss = add2set_loop(ds, H, cp, ccf, iter=i,
                              verbose=verbose,
                              capsize=capsize,
                              currsize=len(s))
        s = s | ss
        i += 1
        if capsize <= len(s): break

    return s

def edgemask(gl,H, cds):
    """given a list of encoded graphs and observed undersampled graph
    H returns a matrix with -1 on diagonal, 0 at the conflicting graph
    combination and encoded graph at non-conflicted
    positions. Furthermore, returns a set of graphs that are in the
    equivalence class of H

    Arguments:
    - `gl`: list of integer encoded graphs
    - `H`: the observed undersampled graph
    """
    n = len(H)
    nl= len(gl)
    s = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)

    for i in xrange(nl):
        for j in xrange(i+1,nl):

            if gl[i] & gl[j]: continue
            if skip_conflict(gl[i], gl[j], cds): continue

            gnum = gl[i] | gl[j]
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
    return mask, s

def ledgemask(gl,H, cds):
    """given a list of encoded graphs and observed undersampled graph
    H returns a matrix with -1 on diagonal, 0 at the conflicting graph
    combination and encoded graph at non-conflicted
    positions. Furthermore, returns a set of graphs that are in the
    equivalence class of H

    Arguments:
    - `gl`: list of integer encoded graphs
    - `H`: the observed undersampled graph
    """
    n = len(H)
    nl= len(gl)
    s = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)

    for i in xrange(nl):
        for j in xrange(i+1,nl):

            if gl[i] & gl[j]: continue
            gnum = gl[i] | gl[j]
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
    return mask, s

def edgeds(mask):
    """construct an edge dictionary from the mask matrix

    Arguments:
    - `mask`:
    """
    ds = {}
    nl = mask.shape[0]
    idx = np.triu_indices(nl,1)
    for i,j in zip(idx[0], idx[1]):
        if mask[i,j]:
            ds[(i,j)] = set()
            conf = set([i,j])
            conf = conf.union(np.where(mask[i,:]==0)[0])
            conf = conf.union(np.where(mask[j,:]==0)[0])
            for k,m in zip(idx[0], idx[1]):
                if not mask[k,m]: continue
                if k in conf: continue
                if m in conf: continue
                if not (k,m) in ds: ds[(i,j)].add(mask[k,m])
            if not ds[(i,j)]: ds.pop((i,j))
    return ds

def edgedsg(mask):
    """construct an edge dictionary from the mask matrix

    Arguments:
    - `mask`:
    """
    ds = {}
    nl = mask.shape[0]
    idx = np.triu_indices(nl,1)
    for i,j in zip(idx[0], idx[1]):
        if mask[i,j]:
            ds[mask[i,j]] = set()
            conf = set([i,j])
            conf = conf.union(np.where(mask[i,:]==0)[0])
            conf = conf.union(np.where(mask[j,:]==0)[0])
            for k,m in zip(idx[0], idx[1]):
                if not mask[k,m]: continue
                if k in conf: continue
                if m in conf: continue
                if not (k,m) in ds: ds[mask[i,j]].add(mask[k,m])
            if not ds[mask[i,j]]: ds.pop(mask[i,j])
    return ds


def quadlister(glist, H, cds):
    n = len(H)
    s = set()
    cache = {}

    def edgemask(gl, H, cds):
        nl= len(gl)
        ss = set()
        mask = np.zeros((nl,nl),'int')
        np.fill_diagonal(mask,-1)
        idx = np.triu_indices(nl,1)
        for i,j in zip(idx[0], idx[1]):
            if gl[i] & gl[j]:
                mask[i,j] = -1
                mask[j,i] = -1
                continue
            if skip_conflict(gl[i], gl[j], cds):
                gnum = gl[i] | gl[j]
                cache[gnum] = False
                continue

            gnum = gl[i] | gl[j]
            if gnum in cache:
                if cache[gnum]:
                    mask[i,j] = gnum
                    mask[j,i] = gnum
            else:
                cache[gnum] = False
                g = bfu.num2CG(gnum,n)
                if not bfu.call_u_conflicts(g, H):
                    if bfu.call_u_equals(g, H): ss.add(gnum)
                    mask[i,j] = gnum
                    mask[j,i] = gnum
                    cache[gnum] = True
        return mask, ss


    def quadmerger(gl, H, cds):
        mask, ss = edgemask(gl, H, cds)
        ds = edgeds(mask)
        #ipdb.set_trace()
        return [[mask[x]]+list(ds[x]) for x in ds], ss

    l = []
    for gl in glist:
        ll, ss = quadmerger(gl, H, cds)
        l.extend(ll)
        s = s|ss

    return l, s


def dceqc(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = confpairs(H)

    glist =  [2**np.arange(n**2)]
    i = 1
    #for i in range(int(np.log2(n**2))):
    while glist != []:
        print i, np.max(map(len,glist)), len(glist)
        glist, ss = quadlister(glist, H, cds)
        s = s|ss
        i += 1
    return s

def quadmerge(gl, H, cds):
    n = len(H)
    l = set()
    s = set()
    mask, ss = edgemask(gl, H, cds)
    s = s | ss
    ds = edgeds(mask)

    #pp = pprint.PrettyPrinter(indent=1)
    #pp.pprint(ds)

    for idx in ds:
        for gn in ds[idx]:
            if mask[idx]&gn: continue
            if skip_conflict(mask[idx], gn, cds): continue
            gnum = mask[idx] | gn
            if gnum in l or gnum in ss: continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)

    return list(l), s

def skip_conflict(g1, g2, ds):
    pss = False
    for ekey in ds:
        if (g1 & ekey) == ekey:
            if ekey in ds and cacheconflicts(g2,ds[ekey]):
                pss = True
                break
    return pss

def edgemask2(gl,H, cds):
    n = len(H)
    nl= len(gl)
    s = set()
    o = set()
    mask = np.zeros((nl,nl),'int')
    np.fill_diagonal(mask,-1)
    for i in xrange(nl):
        for j in xrange(i+1,nl):
            if gl[i] & gl[j]: continue
            gnum = gl[i] | gl[j]
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                if bfu.call_u_equals(g, H): s.add(gnum)
                mask[i,j] = gnum
                mask[j,i] = gnum
            elif bfu.overshoot(g, H): o.add(gnum)
    return mask, s, o # mask, found eqc members, overshoots

def patchmerge(ds, H, cds):
    n = len(H)
    l = set()
    s = set()
    o = set()
    for gkey in ds:
        for num in ds[gkey]:
            if gkey & num: continue
            gnum = gkey | num
            if gnum is s: continue
            if skip_conflictors(gnum, cds): continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)
            elif not gnum in o and bfu.overshoot(g, H): o.add(gnum)
    return l, s, o

def quadmerge2(gl, H, cds):
    n = len(H)

    mask, s, o = edgemask2(gl, H, cds)
    #ipdb.set_trace()
    ds = edgedsg(mask)
    l, ss, oo = patchmerge(ds, H, cds)

    o = o | oo
    s = s | ss

    print 'overshoots: ', len(o)

    return list(l), s

def quadmerge21(gl, H, cds):
    n = len(H)
    l = set()

    mask, ss, o = edgemask2(gl, H, cds)
    idx = np.triu_indices(mask.shape[0], 1)
    print len(o)
    for i in range(len(idx[0])):
        if mask[idx[0][i],idx[1][i]]: l.add(mask[idx[0][i],idx[1][i]])

    return list(l), ss

def dceqclass2(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cp = confpairs(H)
    confs = conflictor_set(H)
    ccf = conflictors(H)

    def prune_loops(gl, H):
        l = []
        for e in gl:
            if e[0] == e[1] and not (e[1] in H[e[0]] and (1,0) in H[e[0]][e[1]]): continue
            l.append(e)
        return l
    edges = gk.edgelist(gk.complement(bfu.num2CG(0,n)))
    edges = prune_loops(edges, H)
    glist = map(lambda x: e2num(x,n),edges)

    #glist =  list(2**np.arange(n**2))
    i = 0
    while glist != []:
        print 2**i, len(glist)
        glist_prev = glist
        glist, ss = quadmerge21(glist, H, confs)
        s = s|ss
        i += 1


    ds = {x: edges for x in glist_prev}

    for j in range(i, len(H)**2):
        ds, ss = add2set_(ds, H, cp, ccf, iter=j, verbose=True)
        s = s | ss
        if not ds: break

    return s


def dceqclass(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = confpairs(H)

    glist =  [0]+list(2**np.arange(n**2))
    i = 1
    while glist != []:
        print i, len(glist)
        glist, ss = quadmerge(glist, H, cds)
        s = s|ss
        i += 1
    return s

def ldceqclass(H,asl=None):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return set()
    n = len(H)
    s = set()
    cds = lconfpairs(H)
    if asl:
        sloops = asl
    else:
        sloops = prune_loops(allsloops(len(H)),H)

    glist =  sloops
    i = 1
    while glist != []:
        print i, len(glist)
        glist, ss = lquadmerge(glist, H, cds)
        s = s|ss
        i += 1
    return s

def lquadmerge(gl, H, cds):
    n = len(H)
    l = set()
    s = set()
    mask, ss = ledgemask(gl, H, cds)
    s = s | ss
    ds = edgeds(mask)

    #pp = pprint.PrettyPrinter(indent=1)
    #pp.pprint(ds)

    for idx in ds:
        for gn in ds[idx]:
            if mask[idx]&gn: continue
            if skip_conflictors(mask[idx], gn, cds): continue
            gnum = mask[idx] | gn
            if gnum in l or gnum in ss: continue
            g = bfu.num2CG(gnum,n)
            if not bfu.call_u_conflicts(g, H):
                l.add(gnum)
                if bfu.call_u_equals(g, H): s.add(gnum)

    return list(l), s

def quadmerge_(glist, H, ds):
    n = len(H)
    gl = set()
    ss = set()
    conflicts = set()
    for gi in combinations(glist, 2):
        if gi[0] & gi[1]: continue
        #if skip_conflict(gi[0], gi[1], ds): continue
        gnum = gi[0] | gi[1]
        if gnum in conflicts: continue
        if skip_conflictors(gnum, ds):
            conflicts.add(gnum)
            continue
        if gnum in gl: continue
        g = bfu.num2CG(gnum,n)
        if not bfu.call_u_conflicts(g, H):
            gl.add(gnum)
            if bfu.call_u_equals(g, H): ss.add(gnum)
        else:
            conflicts.add(gnum)
    return gl, ss

def ecmerge(H):
    """Find all graphs in the same equivalence class with respect to H

    Arguments:
    - `H`: an undersampled graph
    """
    if cmp.isSclique(H):
        print 'not running on superclique'
        return None
    n = len(H)
    s = set()
    ds = confpairs(H)
    ccf = conflictors(H)
    cset = set()
    for e in ccf:
        cset = cset.union(ccf[e])

    glist =  np.r_[[0],2**np.arange(n**2)]
    #glist =  2**np.arange(n**2)

    #glist, ss = quadmerge(glist,H)

    for i in range(int(2*np.log2(n))):
        print i, len(glist)
        glist, ss = quadmerge_(glist,H, cset)
        s = s | ss
    return s

def getrates(g,H):
    n = len(H)
    au = bfu.call_undersamples(g)
    return list(np.where(map(lambda x: x == H, au))[0])

def withrates(s,H):
    n = len(H)
    d = {g:set() for g in s}
    for g in s:
        d[g] = getrates(bfu.num2CG(g,n),H)
    return d

def add2set(gset, elist, H):
    n = len(H)

    s = set()
    ss = set()

    eremove = {e: True for e in elist}

    for gnum in gset:
        g = bfu.num2CG(gnum, n)
        for e in elist:
            if not e[1] in g[e[0]]:
                gk.addanedge(g,e)
                num = bfu.g2num(g)
                if not num in s:
                    au = bfu.call_undersamples(g)
                    if not gk.checkconflict(H, g, au=au):
                        eremove[e] = False
                        s.add(num)
                        if gk.checkequality(H, g, au=au): ss.add(num)
                gk.delanedge(g,e)

    for e in eremove:
        if eremove[e]: elist.remove(e)

    return s, ss, elist

def eqclass_list(H):
    '''
    Find all graphs in the same equivalence class with respect to
    graph H and any undesampling rate.
    '''
    g = {n:{} for n in H}
    s = set()

    edges = gk.edgelist(gk.complement(g))
    #edges = prune_conflicts(H, g, edges)

    gset = set([bfu.g2num(g)])
    for i in range(len(H)**2):
        print i
        gset, ss, edges = add2set(gset, edges, H)
        s = s | ss
        if not edges: break

    return s

def loop2graph(l,n):
    g = {str(i):{} for i in range(1,n+1)}
    for i in range(len(l)-1):
        g[l[i]][l[i+1]] = set([(0,1)])
    g[l[-1]][l[0]] = set([(0,1)])
    return g

def set_loop(loop, graph):
    for i in range(0,len(loop)-1):
        graph[loop[i]][loop[i+1]] = set([(0,1)])
    graph[loop[-1]][loop[0]] = set([(0,1)])

def rotate(l,n): return l[n:] + l[:n]
def get_perm(loop1, loop2, n=None):
    if not n: n = len(loop1)
    basel = [str(i) for i in xrange(1,n+1)]
    diff1 = set(basel) - set(loop1)
    diff2 = set(basel) - set(loop2)
    if loop1[0] in loop2:
        l2 = rotate(loop2, loop2.index(loop1[0]))
    else:
        l2 = loop2
    mp = {}
    for x,y in zip(loop1+list(diff1),l2+list(diff2)):
        mp[x] = y
    return mp

def permute(g, perm):
    gn = {x:{} for x in g}
    for e in g:
        gn[perm[e]] = {perm[x]:g[e][x] for x in g[e]}
    return gn

def permuteAset(s, perm):
    n = len(perm)
    ns = set()
    for e in s:
        ns.add(bfu.g2num(permute(bfu.num2CG(e,n),perm)))
    return ns

def noverlap_loops(loops):
    d = {}
    for l in loops:
        el = []
        for k in loops:
            if not set(l) & set(k):
                el.append(tuple(k))
                #d.setdefault(tuple(l),set()).add(tuple(k))
        d[tuple(l)] = noverlap_loops(el)
    return d

def loop_combinations(loops):
    s = set()
    d = noverlap_loops(loops)
    def dfs_traverse(d, gs):
        if d:
            for e in d:
                dfs_traverse(d[e], gs|set([e]))
        else:
            s.add(frozenset(gs))
    for e in d:
        dfs_traverse(d[e],set([e]))
    return list(s)

def sorted_loops(g):
    l = [x for x in sl.simple_loops(g,0)]
    s = {}
    for e in l:
        s.setdefault(len(e),[]).append(e)
    return s

def loopgroups(g, n=None):
    d = sorted_loops(g)
    if n:
        return loop_combinations(d[n])
    else:
        l=[]
        for key in d:
            l.append(loop_combinations(d[key]))
        return l

def count_loops(n):
    s = 0
    for i in range(1,n+1):
        s += comb(n,i) * math.factorial(i-1)
    return s

def perm_cyclic(l): return [tuple(l[i:]+l[:i]) for i in range(len(l))]
def hashloop(l):
    t = [int(x) for x in l]
    idx = np.argmin(t)
    return tuple(l[idx:]+l[:idx])
def perm_circular_slow2(l):
    s = [tuple(l)]
    c = {}
    c[hashloop(l)] = True
    for e in permutations(l):
        if not hashloop(e) in c:
            s.append(e)
            c[hashloop(e)] = True
    return s

def perm_circular_slow(l):
    s = [tuple(l)]
    c = set(perm_cyclic(l))
    for e in permutations(l):
        if not e in c:
            s.append(e)
            c = c | set(perm_cyclic(e))
    return s
def perm_circular(l, cp=circp):
    r = []
    n = len(l)
    for e in cp[n]:
        r.append([l[i] for i in e])
    return r

def gen_loops(n):
    l = [str(i) for i in range(1,n+1)]
    s = []
    for i in range(1,n+1):
        for e in combinations(l,i):
            s.extend(perm_circular(e))
    return s

def allsloops(n, asl = alloops):
    if asl: return asl[n]
    s = []
    l = gen_loops(n)
    for e in l:
        s.append(bfu.g2num(loop2graph(e,n)))
    return s

def reverse(H, verbose=True, capsize=1000):
    n = len(H)
    s = set()

    g      = gk.superclique(n)
    sloops = set(allsloops(n))

    ds = {bfu.g2num(g): sloops}

    if verbose: print '%3s'%'i'+'%10s'%' graphs'
    i=0
    while ds:
        ds, ss = del_loop(ds, H, iter=i,
                          verbose=verbose,
                          capsize=capsize)
        s = s | ss
        i += 1
        if capsize <= len(s): break

    return s

# ----------------------

def build_loop_step(ds, loop, n, iter=1):
    n2 = n*n +n
    dsr = {}
    s = set()
    ss = set()

    pbar = start_progress_bar(iter, len(ds))

    c = 0

    for gnum in ds:
        c += 1
        pbar.update(c)
        gset = set()
        eset = set()
        for sloop in ds[gnum]:
            num = sloop | gnum
            if not num in s:
                g = bfu.num2CG(num, n)
                s.add(num)
                if bfu.forms_loop(g, loop):
                    ss.add(num)
                else:
                    gset.add((num,sloop))
                    eset.add(sloop)

        for gn,e in gset:
            dsr[gn] = eset - set([e])

    pbar.finish()
    return dsr, ss

def forward_loop_match(loop, n):
    """start with an empty graph and keep adding simple loops until
    the loop is generated at some undersampling rate

    Arguments:
    - `loop`: binary encoding of the loop
    - `n`: number of nodes in the graph
    """
    s = set()
    sloops = allsloops(n)
    ds = {0: sloops}

    i=0
    while ds:
        ds, ss = build_loop_step(ds, loop, n, iter=i)
        s = s | ss
        i += 1

    return s

def delAloop(g, loop):
    n = len(g)
    l = []
    l = [bfu.g2num(ur.loop2graph(s,n)) for s in sls.simple_loops(g,0)]
    l = [num for num in l if not num == loop ]
    print loop, ': ',  l
    return bfu.num2CG(reduce(operator.or_, l),n)

def reverse_loop_match(g,loop):
    """start with a graph and keep removing loops while the loop is still matched

    Arguments:
    - `g`: graph that generates the loop
    - `loop`: the reference loop
    """
    s = set()
    n = len(g)

    def prune(g):
        numh = bfu.g2num(g)
        cannotprune = True
        for l in sls.simple_loops(gk.digonly(g),0):
            gg = delAloop(g,bfu.g2num(loop2graph(l,n)))
            if bfu.forms_loop(gg, loop):
                cannotprune = False
                prune(gg)
        if cannotprune:
            print 'one'
            s.add(g)

    prune(g)
    return s

def reverse_edge_match(g,loop):
    """start with a graph and keep removing loops while the loop is still matched

    Arguments:
    - `g`: graph that generates the loop
    - `loop`: the reference loop
    """
    s = set()
    n = len(g)

    def prune(g):
        numh = bfu.g2num(g)
        cannotprune = True
        for l in gk.edgelist(gk.digonly(g)):
            gk.delanedge(g,l)
            if bfu.forms_loop(g, loop):
                cannotprune = False
                prune(g)
            gk.addanedge(g,l)
        if cannotprune: s.add(bfu.g2num(g))

    prune(g)
    return s

def matchAloop(loop, n):
    """returns a set of minimal graphs that generate this loop

    Arguments:
    - `loop`: binary encoding of the loop
    - `n`: number of nodes in the graph
    """
    s = set()
    l = forward_loop_match(loop,n)
    print len(l)
    for g in l:
        s = s | reverse_edge_match(bfu.num2CG(g,n),loop)

    return s

# ----------------------

def del_loop(ds, H, iter=0, verbose=True, capsize=1000):
    n = len(H)

    dsr = {}
    s = set()
    ss = set()
    print iter,
    for gnum in ds:
        gset = []
        s = set()
        for sloop in ds[gnum]:
            rset = ds[gnum] - set([sloop])
            num = reduce(operator.or_, rset)
            if not num in s:
                g = bfu.num2CG(num, n)
                if bfu.overshoot(g, H):
                    s.add(num)
                    gset.append((num,rset))

        if gset == []:
            print '.',
            ss.add(gnum)

        for gn in gset: dsr[gn[0]] = gn[1]
    print ''
    return dsr, ss

def main():
    g = bfu.ringmore(6,1);
    H = bfu.undersample(g,1);
    ss = liteqclass(H)
    print ss

if __name__ == "__main__":
    main()
