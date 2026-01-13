from networkx import strongly_connected_components
from functools import reduce, wraps
import scipy
import numpy as np
import itertools, copy, time
import random
import sys,os
TOOLSPATH='~/soft/src/dev/craft/gunfolds/tools/'
sys.path.append(os.path.expanduser(TOOLSPATH))

# local packages
import bfutils as bfu
import graphkit as gk
import comparison
import ecj


def isedgesubsetD(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if (0,1) in g2star[n][h]:
                if h in g2[n]:
                    if not (0,1) in g2[n][h]:
                        return False
                else:
                    return False
    return True

def isedgesubset(g2star,g2):
    '''
    check if g2star edges are a subset of those of g2
    '''
    for n in g2star:
        for h in g2star[n]:
            if h in g2[n]:
                #if not (0,1) in g2[n][h]:
                if not g2star[n][h].issubset(g2[n][h]):
                    return False
            else:
                    return False
    return True

def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in range(0, len(l), n):
        yield l[i:i+n]

def purgepath(path, l):
    for i in range(1,len(path)-1):
        l.remove((path[i],path[i+1]))

def next_or_none(it):
    try:
        n = it.next()
    except StopIteration:
        return None
    return n

def try_till_d_path(g,d,gt,order=None):
    k = []
    i = 1
    while not k:
        if order:
            k = [x for x in length_d_paths(g,order(i),d)]
        else:
            k = [x for x in length_d_paths(g,str(i),d)]
        i += 1
        if i > len(g): return []

    ld = []
    for i in range(min(10, len(k))):
        ld.append(len(checkApath(['2']+k[i], gt)))
    idx = np.argmin(ld)
    return k[0]

def try_till_path(g, gt):
    d = len(g)-1
    gx = comparison.graph2nx(g)
    sccl = [x for x in strongly_connected_components(gx)]
    # take the largest
    ln = [len(x) for x in sccl]
    idx = np.argsort(ln)
    d = len(sccl[idx[-1]])-1
    sccl = [sccl[i] for i in idx[::-1]]
    order = [item for sublist in sccl for item in sublist]
    k = []
    while not k:
        if d < 5: return []
        k = try_till_d_path(g,d,gt)
        d -= 1
    return k

def gpurgepath(g,path):
    for i in range(1,len(path)-1):
        del g[path[i]][path[i+1]]

def forks(n,c,el,bl,doempty=lambda x,y: True):
    '''
    INPUT:
    n - single node string
    c - mutable list of children of node n (will be changed as a side effect)
    el - list of edges available for construction (side effect change)
    bl - list of bidirected edges

    OUTPUT:
    l - list of forks
    '''
    l = []
    r = set()
    for p in [x for x in itertools.combinations(c,2)]:
        if doempty(p, bl) and (n,p[0]) in el and (n,p[1]) in el:
            l.append((n,)+p)
            el.remove((n,p[0]))
            el.remove((n,p[1]))
            r.add(p[0])
            r.add(p[1])
    for e in r: c.remove(e)
    return l

def childrenedges(n,c,el,bl):
    '''
    INPUT:
    n - single node string
    c - mutable list of children of node n (will be changed as a side effect)
    el - list of edges available for construction (side effect change)
    bl - list of bidirected edges

    OUTPUT:
    l - list of
    '''
    l = []
    r = set()
    for ch in c:
        if (n,ch) in el:
            l.append((n,ch))
            el.remove((n,ch))
            r.add(ch)
    for e in r: c.remove(e)
    return l

# an empty fork is a for without the bidirected edge
def make_emptyforks(n,c,el,bl):
    return forks(n,c,el,bl,doempty=lambda x,y: not x in y)
def make_fullforks(n,c,el,bl):
    return forks(n,c,el,bl,doempty=lambda x,y: x in y)

def make_longpaths(g, el):
    '''
    el - list of edges that is modified as a side effect
    '''
    l = []
    gc = copy.deepcopy(g)
    for i in range(16):
        k = try_till_path(gc,g)
        if len(k) < 5: break
        if k:
            l.append(('2',)+tuple(k))
            purgepath(l[-1],el)
            gpurgepath(gc,l[-1])
        else:
            break
    return l

def make_allforks_and_rest(g,el,bl, dofullforks=True):
    '''
    el - list of edges that is modified as a side effect
    '''
    l = []
    r = []
    nodes = [n for n in g]
    random.shuffle(nodes)
    for n in nodes:

        c = [e for e in g[n] if (0,1) in g[n][e]]# all children
        if len(c) == 1:
            if (n,c[0]) in el:
                r.append((n,c[0]))
                el.remove((n,c[0]))
        elif len(c) > 1:
            l.extend(make_emptyforks(n,c,el,bl))
            if dofullforks: l.extend(make_fullforks(n,c,el,bl))
            r.extend(childrenedges(n,c,el,bl))
    return l, r

def vedgelist(g, pathtoo=False):
    """ Return a list of tuples for edges of g and forks
    a superugly organically grown function that badly needs refactoring
    """
    l = []
    el = gk.edgelist(g)
    bl = gk.bedgelist(g)

    if pathtoo: l.extend(make_longpaths(g,el))
    l2,r = make_allforks_and_rest(g,el,bl,dofullforks=True)
    l.extend(l2)

    A, singles = makechains(r)

    if singles:
        B, singles = makesinks(singles)
    else:
        B, singles = [], []

    l = longpaths_pick(l)+threedges_pick(l) + A + B + singles
    return l

def twoedges_pick(l):  return [e for e in l if len(e)==2]
def threedges_pick(l): return [e for e in l if len(e)==3]
def longpaths_pick(l): return [e for e in l if len(e)>3 and e[0]=='2']
def makechains(l):
    """ Greedily construct 2 edge chains from edge list
    """
    ends = {e[1]:e for e in l}
    starts = {e[0]:e for e in l}
    r = []
    singles = []
    while l:

        e = l.pop()
        if e[1] in starts and e[0] != e[1] and starts[e[1]] in l:
            r.append(('0', e[0],)+starts[e[1]])
            l.remove(starts[e[1]])
        elif e[0] in ends and e[0] != e[1] and ends[e[0]] in l:
            r.append(('0',)+ends[e[0]]+(e[1],))
            l.remove(ends[e[0]])
        else:
            singles.append(e)
    return r, singles
def makesink(es): return ('1', es[0][0],) + es[1]
def makesinks(l):
    """ Greedily construct 2 edge sinks ( b->a<-c ) from edge list
    """
    sinks = {}
    for e in l:
        if e[1] in sinks:
            sinks[e[1]].append(e)
        else:
            sinks[e[1]] = [e]
    r = []
    singles = []
    for e in sinks:
        if len(sinks[e])>1:
            for es in chunks(sinks[e],2):
                if len(es)==2:
                    r.append(makesink(es))
                else:
                    singles.append(es[0])
        else:
            singles.append(sinks[e][0])
    return r, singles

def selfloop(n,g):
    return n in g[n]

def selfloops(l,g):
    return reduce(lambda x,y: x and y, map(lambda x: selfloop(x,g), l))

def checkbedges(v,bel,g2):
    r = []
    for e in bel:
        if e == tuple(v[1:]) and not selfloops(e, g2):
            r.append(e)
        if e == (v[2],v[1]) and not selfloops(e, g2):
            r.append(e)
    for e in r: bel.remove(e)
    return bel

def checkedge(e, g2):
    if e[0] == e[1]:
        l = [n for n in g2 if n in g2[n]]
        s = set()
        for v in g2[e[0]]: s=s.union(g2[e[0]][v])
        if not (2,0) in s:
            l.remove(e[0])
        return l
    else:
        return g2.keys()

def single_nodes(v,g2):
    """ Returns a list of singleton nodes allowed for merging with v
    """
    l = [(n,n) for n in g2 if not n in v and len(g2[n])>1]
    return l

def checkvedge(v, g2):
    """ Nodes to check to merge the virtual nodes of v ( b<-a->c )
    """
    l = gk.bedgelist(g2)
    if (v[1],v[2]) in l:
        l = single_nodes(v,g2) + checkbedges(v,l,g2)
        for n in v:
            if n in g2[n]: l.append((n,n))
    else:
        l = checkbedges(v,l,g2)
    return list(set(l))

def checkAedge(v, g2):
    """ Nodes to check to merge the virtual nodes of A ( b->a<-c )
    """
    l = []
    # try all pairs but the sources
    for pair in itertools.combinations(g2,2):
        #if pair == (v[1],v[2]): continue
        #if pair == (v[2],v[1]): continue
        l.append(pair)
        l.append(pair[::-1])
    for n in g2:
        l.append((n,n))
    return l

def checkcedge(c, g2):
    """ Nodes to check to merge the virtual nodes of c ( a->b->c )
    """
    l = gk.edgelist(g2)
    return list(set(l))

def checkApath(p, g2):
    sl = [x for x in g2 if selfloop(x,g2)]
    d = len(p) - 2
    l = []
    for n in g2:
        l.extend([tuple(x) for x in length_d_loopy_paths(g2, n, d, p[1:])])
    #k = prunepaths_1D(g2, p, l)
    return l


def isedge(v):  return len(v) == 2 # a->b
def isvedge(v): return len(v) == 3 # b<-a->c
def isCedge(v): return len(v) == 4 and v[0] == '0' # a->b->c
def isAedge(v): return len(v) == 4 and v[0] == '1'# a->c<-b
def isApath(v):  return len(v) >= 4 and v[0] == '2'# a->b->...->z

def checker(n,ee):
    g = bfu.ringmore(n,ee)
    g2 = bfu.increment(g)
    d = checkable(g2)
    t = [len(d[x]) for x in d]
    r = []
    n = len(g2)
    ee= len(gk.edgelist(g2))
    for i in range(1,len(t)):
        r.append(sum(np.log10(t[:i])) - ee*np.log10(n))
    return r

def checkerDS(n,ee):
    g = bfu.ringmore(n,ee)
    g2 = bfu.increment(g)
    gg = checkable(g2)
    d,p,idx = conformanceDS(g2,gg,gg.keys())
    t = [len(x) for x in p]
    r = []
    n = len(g2)
    ee= len(gk.edgelist(g2))
    for i in range(1,len(t)):
        r.append(sum(np.log10(t[:i])) - ee*np.log10(n))
    return r

def fordens(n,denslist, repeats=100):
    rl={}
    for d in denslist:
        ee = bfu.dens2edgenum(d,n)
        l=[checker(n,ee)[-1] for i in range(repeats)]
        rl[d] = (round(scipy.mean(l),3),round(scipy.std(l),3))
    return rl

def fordensDS(n,denslist, repeats=100):
    rl={}
    for d in denslist:
        print(d)
        ee = bfu.dens2edgenum(d,n)
        l=[checkerDS(n,ee)[-1] for i in range(repeats)]
        rl[d] = (round(scipy.mean(l),3),round(scipy.std(l),3))
    return rl

def checkable(g2):
    d = {}
    g = cloneempty(g2)
    vlist = vedgelist(g2,pathtoo=False)
    for v in vlist:
        if isvedge(v):
            d[v] = checkvedge(v,g2)
        elif isCedge(v):
            d[v] = checkcedge(v,g2)
        elif isAedge(v):
            d[v] = checkAedge(v,g2)
        elif isApath(v):
            d[v] = checkApath(v,g2)
        else:
            d[v] = checkedge(v,g2)

    # # check if some of the otherwise permissible nodes still fail
    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]

    for e in d:
        adder, remover = f[edge_function_idx(e)]
        checks_ok = c[edge_function_idx(e)]
        for n in d[e]:
            if not checks_ok(e,n,g,g2):
                d[e].remove(n)

    return d

def inorder_check2(e1, e2, j1, j2, g2, f=[], c=[]):
    g = cloneempty(g2) # the graph to be used for checking

    if f==[]:
        f = [(add2edges,del2edges,mask2edges),
             (addavedge,delavedge,maskavedge),
             (addacedge,delacedge,maskaCedge),
             (addaAedge,delaAedge,maskaAedge),
             (addapath,delapath,maskapath)]

    if c==[]:
        c = [ok2add2edges,
             ok2addavedge,
             ok2addacedge,
             ok2addaAedge,
             ok2addapath]

    adder, remover, masker = f[edge_function_idx(e1)]
    checks_ok = c[edge_function_idx(e2)]

    d = {}
    s1 = set()
    s2 = set()
    for c1 in j1: # for each connector
        mask = adder(g,e1,c1)
        d[c1] = set()
        for c2 in j2:
            if checks_ok(e2,c2,g,g2):
               d[c1].add(c2)
               s2.add(c2)
        remover(g,e1,c1,mask)
        if d[c1]: s1.add(c1)
    return d,s1,s2

def check3(e1, e2, e3, j1, j2, j3, g2, f=[], c=[]):
    g = cloneempty(g2) # the graph to be used for checking
    if f==[]:
        f = [(add2edges,del2edges,mask2edges),
             (addavedge,delavedge,maskavedge),
             (addacedge,delacedge,maskaCedge),
             (addaAedge,delaAedge,maskaAedge),
             (addapath,delapath,maskapath)]
    if c==[]:
        c = [ok2add2edges,
             ok2addavedge,
             ok2addacedge,
             ok2addaAedge,
             ok2addapath]

    adder1, remover1, masker1 = f[edge_function_idx(e1)]
    adder2, remover2, masker2 = f[edge_function_idx(e2)]

    checks_ok2 = c[edge_function_idx(e2)]
    checks_ok3 = c[edge_function_idx(e3)]

    d = {}
    s1 = set()
    s2 = set()
    s3 = set()

    for c1 in j1: # for each connector
        mask1 = adder1(g,e1,c1)
        append_set1 = False
        for c2 in j2:
            append_set2 = False
            if checks_ok2(e2,c2,g,g2):
                mask2 = adder2(g,e2,c2)
                for c3 in j3:
                    if checks_ok3(e3,c3,g,g2):
                        append_set1 = append_set2 = True
                        s3.add(c3)
                remover2(g,e2,c2,mask2)
            if append_set2: s2.add(c2)
        if append_set1: s1.add(c1)
        remover1(g,e1,c1,mask1)
    return s1, s2, s3

def del_empty(d):
    l = [e for e in d]
    for e in l:
        if d[e]==set(): del d[e]
    return d
def inorder_checks(g2, gg):
    #idx = np.argsort([len(gg[x]) for x in gg.keys()])
    #ee = [gg.keys()[i] for i in idx] # to preserve the order
    ee = [e for e in gg] # to preserve the order
    #cds = conformanceDS(g2, ee)
    #oo = new_order(g2, ee, repeats=100, cds=None)
    #ee = oo[0]
    random.shuffle(ee)
    d = {} # new datastructure
    d[ee[0]] = {('0'):gg[ee[0]]}
    for i in range(len(ee)-1):
        d[ee[i+1]] = del_empty(inorder_check2(ee[i], ee[i+1],
                                              gg[ee[i]], gg[ee[i+1]], g2)[0])
    return ee, d

def cloneempty(g): return {n:{} for n in g} # return a graph with no edges

def ok2addanedge1(s, e, g, g2,rate=1):
    """
    s - start,
    e - end
    """
    # directed edges
    # self-loop
    if s == e and not e in g2[s]: return False
    for u in g: # Pa(s) -> e
        if s in g[u] and not (e in g2[u] and (0,1) in g2[u][e]):
            return False
    for u in g[e]: # s -> Ch(e)
        if not (u in g2[s] and (0,1) in g2[s][u]):return False
    # bidirected edges
    for u in g[s]: # e <-> Ch(s)
        if u!=e and not (u in g2[e] and (2,0) in g2[e][u]):return False
    return True

def ok2addanedge2(s, e, g, g2, rate=1):
    mask = addanedge(g,(s,e))
    value = bfu.undersample(g,rate) == g2
    delanedge(g,(s,e),mask)
    return value

def ok2addanedge_sub(s, e, g, g2, rate=1):
    mask = addanedge(g,(s,e))
    value = isedgesubset(bfu.undersample(g,rate),g2)
    delanedge(g,(s,e),mask)
    return value

def ok2addanedge(s, e, g, g2, rate=1):
    f = [ok2addanedge1, ok2addanedge2]
    return f[min([1,rate-1])](s,e,g,g2,rate=rate)

def ok2addanedge_(s, e, g, g2, rate=1):
    f = [ok2addanedge1, ok2addanedge_sub]
    return f[min([1,rate-1])](s,e,g,g2,rate=rate)


def ok2add2edges(e,p,g,g2): return edge_increment_ok(e[0],p,e[1],g,g2)

def maskanedge(g,e): return [e[1] in g[e[0]]]
def mask2edges(g,e,p): return [p in g[e[0]], e[1] in g[p]]
def maskavedge(g,e,p):
    return [p[0] in g[e[0]], p[1] in g[e[0]],
            e[1] in g[p[0]], e[2] in g[p[1]]]
def maskaAedge(g,e,p):
    return [p[0] in g[e[1]], p[1] in g[e[2]],
            e[3] in g[p[0]], e[3] in g[p[1]]]
def maskaCedge(g,e,p):
    return [p[0] in g[e[1]], e[2] in g[p[0]],
            p[1] in g[e[2]], e[3] in g[p[1]]]
def maskapath(g,e,p):
    mask = []
    for i in range(len(p)):
        mask.append(p[i] in g[e[i+1]])
        mask.append(e[i+2] in g[p[i]])
    return mask
def maskaVpath(g,e,p):
    mask = []
    mask.extend([p[0] in g[e[0]], e[1] in g[p[-1]]])
    for i in range(1,len(p)):
        mask.append(p[i] in g[p[i-1]])
    return mask

def addanedge(g,e):
    '''
    add edge e[0] -> e[1] to g
    '''
    mask = maskanedge(g,e)
    g[e[0]][e[1]] =  set([(0,1)])
    return mask
def delanedge(g,e,mask):
    '''
    delete edge e[0] -> e[1] from g if it was not there before
    '''
    if not mask[0]: g[e[0]].pop(e[1], None)


def add2edges(g,e,p):
    '''
    break edge e[0] -> e[1] into two pieces
    e[0] -> p and p -> e[1]
    and add them to g
    '''
    mask = mask2edges(g,e,p)
    g[e[0]][p] = g[p][e[1]] = set([(0,1)])
    return mask

def del2edges(g,e,p,mask):
    '''
    restore the graph as it was before adding e[0]->p and p->e[1]
    '''
    if not mask[0]: g[e[0]].pop(p, None)
    if not mask[1]: g[p].pop(e[1], None)

def ok2addavedge(e,p,g,g2):
    if p[1] == e[0]:
        if p[0] != p[1] and p[0] != e[2] and not (e[2] in g2[p[0]] and (2,0) in g2[p[0]][e[2]]):
                return False
        if p[0] == p[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
        if p[0] == e[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if p[0] == e[0]:
        if p[0] != p[1] and p[1] != e[1] and not (e[1] in g2[p[1]] and (2,0) in g2[p[1]][e[1]]):
                return False
        if p[0] == p[1] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
        if p[1] == e[2] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if p[0] == e[1] and p[1] == e[2] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False
    if p[0] == e[2] and not (e[1] in g2[p[1]] and (0,1) in g2[p[1]][e[1]]):
            return False
    if p[1] == e[1] and not (e[2] in g2[p[0]] and (0,1) in g2[p[0]][e[2]]):
            return False
    if p[0] == p[1] == e[0] and not (e[2] in g2[e[1]] and (2,0) in g2[e[1]][e[2]]):
            return False

    if not edge_increment_ok(e[0],p[0],e[1],g,g2): return False
    if not edge_increment_ok(e[0],p[1],e[2],g,g2): return False

    return  True

def addavedge(g,v,b):
    mask = maskavedge(g,v,b)
    g[v[0]][b[0]] = g[v[0]][b[1]] = g[b[0]][v[1]] = g[b[1]][v[2]] = set([(0,1)])
    return mask

def delavedge(g,v,b,mask):
    if not mask[0]: g[v[0]].pop(b[0], None)
    if not mask[1]: g[v[0]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[1], None)
    if not mask[3]: g[b[1]].pop(v[2], None)

def ok2addaAedge(e,p,g,g2):
    if p[1] == e[1] and not (p[0] in g2[e[2]] and (0,1) in g2[e[2]][p[0]]): return False
    if p[0] == e[2] and not (p[1] in g2[e[1]] and (0,1) in g2[e[1]][p[1]]): return False

    if not edge_increment_ok(e[1],p[0],e[3],g,g2): return False
    if not edge_increment_ok(e[2],p[1],e[3],g,g2): return False

    return True

def addaAedge(g,v,b):
    mask = maskaAedge(g,v,b)
    g[v[1]][b[0]] = g[v[2]][b[1]] = g[b[0]][v[3]] = g[b[1]][v[3]] = set([(0,1)])
    return mask

def delaAedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[v[2]].pop(b[1], None)
    if not mask[2]: g[b[0]].pop(v[3], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def cleanedges(e,p,g, mask):
    i = 0
    for m in mask:
        if not m[0]: g[e[i+1]].pop(p[i], None)
        if not m[1]: g[p[i]].pop(e[i+2], None)
        i += 1
def cleanVedges(g, e,p, mask):

    if mask:
        if not mask[0]: g[e[0]].pop(p[0], None)
        if not mask[1]: g[p[-1]].pop(e[1], None)

        i = 0
        for m in mask[2:]:
            if not m: g[p[i]].pop(p[i+1], None)
            i += 1

def ok2addapath(e,p,g,g2):
    mask = []
    for i in range(len(p)):
        if not edge_increment_ok(e[i+1],p[i],e[i+2],g,g2):
            cleanedges(e,p,g,mask)
            return False
        mask.append(add2edges(g,(e[i+1],e[i+2]),p[i]))
    cleanedges(e,p,g,mask)
    return True

def ok2addaVpath(e,p,g,g2,rate=2):
    mask = addaVpath(g,e,p)
    if not isedgesubset(bfu.undersample(g,rate), g2):
        cleanVedges(g,e,p,mask)
        return False
    cleanVedges(g,e,p,mask)
    return True

def ok2addapath1(e,p,g,g2):
    for i in range(len(p)):
        if not edge_increment_ok(e[i+1],p[i],e[i+2],g,g2):
            return False
    return True

def addapath(g,v,b):

    mask = maskapath(g,v,b)
    s = set([(0,1)])
    for i in range(len(b)):
        g[v[i+1]][b[i]] = g[b[i]][v[i+2]] = s

    return mask

def addaVpath(g,v,b):
    mask = maskaVpath(g,v,b)
    s = set([(0,1)])
    l = [v[0]] + list(b) + [v[1]]
    for i in range(len(l)-1):
        g[l[i]][l[i+1]] = s
    return mask

def delaVpath(g, v, b, mask):
    cleanVedges(g, v, b, mask)

def delapath(g, v, b, mask):
    for i in range(len(b)):
        if not mask[2*i]: g[v[i+1]].pop(b[i], None)
        if not mask[2*i+1]:g[b[i]].pop(v[i+2], None)

def prunepaths_1D(g2, path, conn):
    c = []
    g = cloneempty(g2)
    for p in conn:
        mask = addapath(g,path,p)
        if isedgesubset(bfu.increment(g), g2): c.append(tuple(p))
        delapath(g,path,p,mask)
    return c

def ok2addacedge(e,p,g,g2):

    if p[0] == p[1]:
        if not e[2] in g2[e[2]]: return False
        if not p[0] in g2[p[0]]: return False
        if not (e[3] in g2[e[1]] and (0,1) in g2[e[1]][e[3]]): return False

    if not edge_increment_ok(e[1],p[0],e[2],g,g2): return False
    if not edge_increment_ok(e[2],p[1],e[3],g,g2): return False

    return True

def addacedge(g,v,b): # chain
    mask = maskaCedge(g,v,b)
    g[v[1]][b[0]] = g[v[2]][b[1]] = g[b[0]][v[2]] = g[b[1]][v[3]] = set([(0,1)])
    return mask

def delacedge(g,v,b,mask):
    if not mask[0]: g[v[1]].pop(b[0], None)
    if not mask[1]: g[b[0]].pop(v[2], None)
    if not mask[2]: g[v[2]].pop(b[1], None)
    if not mask[3]: g[b[1]].pop(v[3], None)

def rotate(l): return l[1:] + l[:1] # rotate a list
def density(g): return len(gk.edgelist(g))/np.double(len(g)**2)
def udensity(g): return (len(gk.edgelist(g))+len(gk.bedgelist(g))/2.)/np.double(len(g)**2 + len(g)*(len(g)-1)/2.)

def esig(l,n):
    '''
    turns edge list into a hash string
    '''
    z = len(str(n))
    n = map(lambda x: ''.join(map(lambda y: y.zfill(z),x)), l)
    n.sort()
    n = ''.join(n[::-1])
    return int('1'+n)

def gsig(g):
    '''
    turns input graph g into a hash string using edges
    '''
    return bfu.g2num(g)

def signature(g, edges): return (gsig(g),esig(edges,len(g)))

def memo(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def memo1(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = gsig(args[0])             # Signature: g
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return cache[s]               # Return the cached solution
    return wrap

def memo2(func):
    cache = {}                        # Stored subproblem solutions
    @wraps(func)                      # Make wrap look like func
    def wrap(*args):                  # The memoized wrapper
        s = signature(args[0],args[2])# Signature: g and edges
        if s not in cache:            # Not already computed?
            cache[s] = func(*args)    # Compute & cache the solution
        return set()#cache[s]               # Return the cached solution
    return wrap

def eqsearch(g2, rate=1):
    '''Find  all  g  that are also in  the equivalence
    class with respect to g2 and the rate.
    '''

    s = set()
    noop = set()

    @memo1
    def addnodes(g,g2,edges):
        if edges:
            masks  = []
            for e in edges:
                if ok2addanedge_(e[0],e[1],g,g2,rate=rate):
                    masks.append(True)
                else:
                    masks.append(False)
            nedges = [edges[i] for i in range(len(edges)) if masks[i]]
            n = len(nedges)
            if n:
                for i in range(n):
                    mask = addanedge(g,nedges[i])
                    if bfu.undersample(g,rate) == g2: s.add(bfu.g2num(g))
                    addnodes(g,g2,nedges[:i]+nedges[i+1:])
                    delanedge(g,nedges[i],mask)
                return s
            else:
                return noop
        else:
            return noop


    g = cloneempty(g2)
    edges = gk.edgelist(gk.complement(g))
    addnodes(g,g2,edges)
    return s


def supergraphs_in_eq(g, g2, rate=1):
    '''Find  all supergraphs of g  that are also in  the same equivalence
    class with respect to g2 and the rate.
    Currently works only for bfu.undersample by 1
    '''
    if bfu.undersample(g,rate) != g2:
        raise ValueError('g is not in equivalence class of g2')

    s = set()

    def addnodes(g,g2,edges):
        if edges:
            masks  = []
            for e in edges:
                if ok2addanedge(e[0],e[1],g,g2,rate=rate):
                    masks.append(True)
                else:
                    masks.append(False)
            nedges = [edges[i] for i in range(len(edges)) if masks[i]]
            n = len(nedges)
            if n:
                for i in range(n):
                    mask = addanedge(g,nedges[i])
                    s.add(bfu.g2num(g))
                    addnodes(g,g2,nedges[:i]+nedges[i+1:])
                    delanedge(g,nedges[i],mask)

    edges = gk.edgelist(gk.complement(g))
    addnodes(g,g2,edges)
    return s

def edge_backtrack2g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print('Superclique - any SCC with GCD = 1 fits')
        return set([-1])

    single_cache = {}

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            ln = [n for n in g2]
            for n in ln:
                if (n,e) in single_cache: continue
                mask = add2edges(g,e,n)
                if isedgesubset(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and bfu.increment(r)==g2:
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                del2edges(g,e,n,mask)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:
            mask = add2edges(g,e,n)
            if not isedgesubset(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def edge_backtrack2g1_directed(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print('Superclique - any SCC with GCD = 1 fits')
        return set([-1])

    single_cache = {}

    def edgeset(g):
        return set(gk.edgelist(g))
    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            e = edges.pop()
            ln = [n for n in g2]
            for n in ln:
                if (n,e) in single_cache: continue
                mask = add2edges(g,e,n)
                if isedgesubsetD(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and edgeset(bfu.increment(r))==edgeset(g2):
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements in eqclass')
                del2edges(g,e,n,mask)
            edges.append(e)
        else:
            return g
    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:
            mask = add2edges(g,e,n)
            if not isedgesubsetD(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def g22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print('Superclique - any SCC with GCD = 1 fits')
        return set([-1])

    single_cache = {}

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            if bfu.increment(g) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                return g
            e = edges[0]
            for n in g2:

                if (n,e) in single_cache: continue
                if not edge_increment_ok(e[0],n,e[1],g,g2): continue

                mask = add2edges(g,e,n)
                r = nodesearch(g,g2,edges[1:],s)
                del2edges(g,e,n,mask)

        elif bfu.increment(g)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements in eqclass')
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in g2:

            mask = add2edges(g,e,n)
            if not isedgesubset(bfu.increment(g), g2):
                single_cache[(n,e)] = False
            del2edges(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def backtrack_more(g2, rate=1, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print('Superclique - any SCC with GCD = 1 fits')
        return set([-1])

    single_cache = {}
    if rate == 1:
        ln = [n for n in g2]
    else:
        ln = []
        for x in itertools.combinations_with_replacement(g2.keys(),rate):
            ln.extend(itertools.permutations(x,rate))
        ln = set(ln)

    @memo # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            if bfu.undersample(g,rate) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                return g
            e = edges[0]
            for n in ln:

                if (n,e) in single_cache: continue
                if not ok2addaVpath(e, n, g, g2, rate=rate): continue

                mask = addaVpath(g,e,n)
                r = nodesearch(g,g2,edges[1:],s)
                delaVpath(g,e,n,mask)

        elif bfu.undersample(g,rate)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements in eqclass')
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    edges = gk.edgelist(g2)
    random.shuffle(edges)
    g = cloneempty(g2)

    for e in edges:
        for n in ln:

            mask = addaVpath(g,e,n)
            if not isedgesubset(bfu.undersample(g,rate), g2):
                single_cache[(n,e)] = False
            delaVpath(g,e,n,mask)

    s = set()
    try:
        nodesearch(g,g2,edges,s)
    except ValueError:
        s.add(0)
    return s

def backtrackup2u(H,umax=2):
    s = set()
    for i in range(1,umax+1):
        s = s | backtrack_more(H,rate=i)
    return s

def vg22g1(g2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print('Superclique - any SCC with GCD = 1 fits')
        return set([-1])

    f = [(add2edges, del2edges),
         (addavedge,delavedge),
         (addacedge,delacedge),
         (addaAedge,delaAedge),
         (addapath,delapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]
    @memo2 # memoize the search
    def nodesearch(g, g2, edges, s):
        if edges:
            #key, checklist = edges.popitem()
            key = random.choice(edges.keys())
            checklist = edges.pop(key)
            adder, remover = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]
            for n in checklist:
                mask = adder(g,key,n)
                if isedgesubset(bfu.increment(g), g2):
                    r = nodesearch(g,g2,edges,s)
                    if r and bfu.increment(r)==g2:
                        s.add(bfu.g2num(r))
                        if capsize and len(s)>capsize:
                            raise ValueError('Too many elements')
                remover(g,key,n,mask)
            edges[key] = checklist
        else:
            return g

    # find all directed g1's not conflicting with g2
    n = len(g2)
    chlist = checkable(g2)
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g,g2,chlist,s)
    except ValueError:
        s.add(0)
    return s

def edge_function_idx(edge):
    return min(4,len(edge))-2+min(max(3,len(edge))-3,1)*int(edge[0])

def v2g22g1(g2, capsize=None, verbose=True):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print('Superclique - any SCC with GCD = 1 fits')
        return set([-1])

    f = [(add2edges,del2edges,mask2edges),
         (addavedge,delavedge,maskavedge),
         (addacedge,delacedge,maskaCedge),
         (addaAedge,delaAedge,maskaAedge),
         (addapath,delapath,maskapath)]
    c = [ok2add2edges,
         ok2addavedge,
         ok2addacedge,
         ok2addaAedge,
         ok2addapath]

    def predictive_check(g,g2,pool,checks_ok, key):
        s = set()
        for u in pool:
            if not checks_ok(key,u,g,g2): continue
            s.add(u)
        return s

    @memo2 # memoize the search
    def nodesearch(g, g2, order, inlist, s, cds, pool, pc):
        if order:
            if bfu.increment(g) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                s.update(supergraphs_in_eq(g, g2))
                return g

            key = order[0]
            if pc:
                tocheck = [x for x in pc if x in cds[len(inlist)-1][inlist[0]]]
            else:
                tocheck = cds[len(inlist)-1][inlist[0]]

            if len(order) > 1:
                kk = order[1]
                pc = predictive_check(g,g2,pool[len(inlist)],
                                      c[edge_function_idx(kk)],kk)
            else:
                pc = set()

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            for n in tocheck:
                if not checks_ok(key,n,g,g2): continue
                masked = np.prod(masker(g,key,n))
                if masked:
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                else:
                    mask = adder(g,key,n)
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                    remover(g,key,n,mask)

        elif bfu.increment(g)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements')
            return g

    @memo2 # memoize the search
    def nodesearch0(g, g2, order, inlist, s, cds):

        if order:
            key = order.pop(0)
            tocheck = cds[len(inlist)-1][inlist[0]]

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            if len(tocheck) > 1:
                for n in tocheck:
                    if not checks_ok(key,n,g,g2): continue
                    mask = masker(g,key,n)
                    if not np.prod(mask):
                        mask = adder(g,key,n)
                        r = nodesearch0(g,g2,order, [n]+inlist, s, cds)
                        if r and bfu.increment(r)==g2:
                            s.add(bfu.g2num(r))
                            if capsize and len(s)>capsize:
                                raise ValueError('Too many elements')
                        remover(g,key,n,mask)
                    else:
                        r = nodesearch0(g,g2,order, [n]+inlist, s, cds)
                        if r and bfu.increment(r)==g2:
                            s.add(bfu.g2num(r))
                            if capsize and len(s)>capsize:
                                raise ValueError('Too many elements')
            elif tocheck:
                (n,) = tocheck
                mask = adder(g,key,n)
                r = nodesearch0(g,g2, order, [n]+inlist, s, cds)
                if r and bfu.increment(r) == g2:
                    s.add(bfu.g2num(r))
                    if capsize and len(s)>capsize:
                        raise ValueError('Too many elements')
                remover(g,key,n,mask)

            order.insert(0,key)

        else:
            return g

    # find all directed g1's not conflicting with g2

    startTime = int(round(time.time() * 1000))
    gg = checkable(g2)

    idx = np.argsort([len(gg[x]) for x in gg.keys()])
    keys = [gg.keys()[i] for i in idx]

    cds, order, idx = conformanceDS(g2, gg, keys)
    endTime = int(round(time.time() * 1000))
    if verbose:
        print("precomputed in {:10} seconds".format(round((endTime-startTime)/1000.,3)))
    if 0 in [len(x) for x in order]:
        return set()
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g, g2, [keys[i] for i in idx], ['0'], s, cds, order, set())
        #nodesearch0(g, g2, [gg.keys()[i] for i in idx], ['0'], s, cds)
    except ValueError as e:
        print(e)
        s.add(0)
    return s

def backtrack_more2(g2, rate=2, capsize=None):
    '''
    computes all g1 that are in the equivalence class for g2
    '''
    if ecj.isSclique(g2):
        print('Superclique - any SCC with GCD = 1 fits')
        return set([-1])

    f = [(addaVpath,delaVpath,maskaVpath)]
    c = [ok2addaVpath]

    def predictive_check(g,g2,pool,checks_ok, key):
        s = set()
        for u in pool:
            if not checks_ok(key,u,g,g2,rate=rate): continue
            s.add(u)
        return s

    @memo2 # memoize the search
    def nodesearch(g, g2, order, inlist, s, cds, pool, pc):
        if order:
            if bfu.undersample(g,rate) == g2:
                s.add(bfu.g2num(g))
                if capsize and len(s)>capsize:
                    raise ValueError('Too many elements')
                s.update(supergraphs_in_eq(g, g2, rate=rate))
                return g

            key = order[0]
            if pc:
                tocheck = [x for x in pc if x in cds[len(inlist)-1][inlist[0]]]
            else:
                tocheck = cds[len(inlist)-1][inlist[0]]

            if len(order) > 1:
                kk = order[1]
                pc = predictive_check(g,g2,pool[len(inlist)],
                                      c[edge_function_idx(kk)],kk)
            else:
                pc = set()

            adder, remover, masker = f[edge_function_idx(key)]
            checks_ok = c[edge_function_idx(key)]

            for n in tocheck:
                if not checks_ok(key,n,g,g2,rate=rate): continue
                masked = np.prod(masker(g,key,n))
                if masked:
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                else:
                    mask = adder(g,key,n)
                    nodesearch(g,g2,order[1:], [n]+inlist, s, cds, pool, pc)
                    remover(g,key,n,mask)

        elif bfu.undersample(g,rate)==g2:
            s.add(bfu.g2num(g))
            if capsize and len(s)>capsize:
                raise ValueError('Too many elements')
            return g

    # find all directed g1's not conflicting with g2

    startTime = int(round(time.time() * 1000))
    ln = [x for x in itertools.permutations(g2.keys(),rate)] + \
         [(n,n) for n in g2]
    gg = {x:ln for x in gk.edgelist(g2)}
    keys = gg.keys()
    cds, order, idx = conformanceDS(g2, gg, gg.keys(), f=f, c=c)
    endTime = int(round(time.time() * 1000))
    print("precomputed in {:10} seconds".format(round((endTime-startTime)/1000.,3)))
    if 0 in [len(x) for x in order]:
        return set()
    g = cloneempty(g2)

    s = set()
    try:
        nodesearch(g, g2, [keys[i] for i in idx], ['0'], s, cds, order, set())
    except ValueError as e:
        print(e)
        s.add(0)
    return s


def unionpool(idx, cds):
    s = set()
    for u in cds[idx]:
        for v in cds[idx][u]:
            s = s.union(cds[idx][u][v])
    return s

def conformanceDS(g2, gg, order, f=[], c=[]):
    CDS = {}
    pool = {}

    CDS[0] = set(gg[order[0]])
    pool = [set(gg[order[i]]) for i in range(len(order))]

    for x in itertools.combinations(range(len(order)),2):

        d, s_i1, s_i2 = inorder_check2(order[x[0]], order[x[1]],
                                       pool[x[0]], pool[x[1]],
                                       g2, f=f, c=c)

        pool[x[0]] = pool[x[0]].intersection(s_i1)
        pool[x[1]] = pool[x[1]].intersection(s_i2)

        d = del_empty(d)
        if not x[1] in CDS:
            CDS[x[1]] = {}
            CDS[x[1]][x[0]] = d
        else:
            CDS[x[1]][x[0]] = d
    if density(g2) > 0.35:
        itr3 = [x for x in itertools.combinations(range(len(order)),3)]
        for x in random.sample(itr3, min(10,np.int(scipy.misc.comb(len(order),3)))):
            s1, s2, s3 = check3(order[x[0]], order[x[1]], order[x[2]],
                                pool[x[0]], pool[x[1]], pool[x[2]],
                                g2, f=f, c=c)

            pool[x[0]] = pool[x[0]].intersection(s1)
            pool[x[1]] = pool[x[1]].intersection(s2)
            pool[x[2]] = pool[x[2]].intersection(s3)

    return prune_sort_CDS(CDS, pool)

def prune_modify_CDS(cds, pool):
    ds = {}
    ds[0]={}
    ds[0]['0'] = pool[0]
    for i in range(1,len(pool)):
        ds[i] = {}
        for j in cds[i].keys():
            for e in pool[i-1].intersection(cds[i][j].keys()):
                ds[i][e] = pool[i].intersection(cds[i][j][e])
    return ds, pool, range(len(pool))

def prune_sort_CDS(cds, pool):

    idx = np.argsort([len(x) for x in pool])
    p = [pool[i] for i in idx]

    ds = {}
    ds[0]={}
    ds[0]['0'] = pool[idx[0]]

    for i in range(1,len(idx)):
        ds[i] = {}
        for j in range(i):
            if idx[j] > idx[i]:
                dd = invertCDSelement(cds[idx[j]][idx[i]])
            else:
                dd = cds[idx[i]][idx[j]]
            for e in pool[idx[i-1]].intersection(dd.keys()):
                ds[i][e] = pool[idx[i]].intersection(dd[e])

    return ds, p, idx

def invertCDSelement(d_i):
    d = {}
    for e in d_i:
        for v in d_i[e]:
            if v in d:
                d[v].add(e)
            else:
                d[v]=set([e])
    return d

def conformant(cds, inlist):

    if inlist[len(inlist)-2] in cds[len(inlist)-1][0]:
        s = cds[len(inlist)-1][0][inlist[len(inlist)-2]]
    else:
        return set()
    for i in range(1,len(inlist)-1):
        if inlist[len(inlist)-i-2] in cds[len(inlist)-1][i]:
            s = s.intersection(cds[len(inlist)-1][i][inlist[len(inlist)-i-2]])
        else:
            return set()
    return s

def length_d_paths(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    def recurse(G, s, d, path=[]):
        if d == 0:
            yield path
            return

        for u in G[s]:
            if G[s][u] == set([(2,0)]) or u in path: continue
            for v in recurse(G, u, d-1, path+[u]):
                yield v

    for u in recurse(G, s, d, [s]):
            yield u

def edge_increment_ok(s,m,e,g,g2):
    """
    s - start,
    m - middle,
    e - end
    """
    # bidirected edges
    for u in g[s]:
        if u!=m and not (m in g2[u] and (2,0) in g2[u][m]):return False

    # directed edges
    if s == e:
        if not (m in g2[m] and (0,1) in g2[m][m]):return False
        if not (s in g2[s] and (0,1) in g2[s][s]):return False
    for u in g[m]:
        if not (u in g2[s] and (0,1) in g2[s][u]):return False
        # bidirected edges
        if u!=e and not (e in g2[u] and (2,0) in g2[u][e]):return False
    for u in g[e]:
        if not (u in g2[m] and (0,1) in g2[m][u]):return False

    for u in g:
        if s in g[u] and not (m in g2[u] and (0,1) in g2[u][m]):
            return False
        if m in g[u] and not (e in g2[u] and (0,1) in g2[u][e]):
            return False

    return True


def length_d_loopy_paths(G, s, dt, p):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    g1 = cloneempty(G)

    def recurse(g, g2, s, d, path=[]):

        if edge_increment_ok(p[-d-2],s,p[-d-1],g,g2):

            if d == 0:
                yield path
                return

            mask = add2edges(g,(p[-d-2],p[-d-1]),s)
            for u in g2[s]:
                if g2[s][u] == set([(2,0)]): continue
                for v in recurse(g, g2, u, d-1, path+[u]):
                    yield v
            del2edges(g,(p[-d-2],p[-d-1]),s, mask)

    for u in recurse(g1, G, s, dt-1, [s]):
            yield u
