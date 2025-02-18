import scipy
from itertools import combinations
from testgraphs import *

def walk(G, s, S=set()):
    P, Q = dict(), set()
    P[s] = None
    Q.add(s)
    while Q:
        u = Q.pop()
        for v in G[u].difference(P,S):
            Q.add(v)
            P[v] = u
    return P

def traverse(G, s, qtype=set):
    S, Q = set(), qtype()
    Q.add(s)
    while Q:
        u = Q.pop()
        if u in S: continue
        S.add(u)
        for v in G[u]:
            Q.add(v)
        yield u

def dfs_topsort(G):
    S, res = set(), []
    def recurse(u):
        if u in S: return
        S.add(u)
        for v in G[u]:
            recurse(v)
        res.append(u)
    for u in G:
        recurse(u)
    res.reverse()
    return res

def tr(G):                      # Transpose (rev. edges of) G
    GT = {}
    for u in G: GT[u] = set()   # Get all the nodes in there
    for u in G:
        for v in G[u]:
            GT[v].add(u)        # Add all reverse edges
    return GT

def scc(G):                   # Kosaraju's algorithm
    GT = tr(G)                # Get the transposed graph
    sccs, seen = [], set()
    for u in dfs_topsort(G):   # DFS starting points
        if u in seen: continue # Ignore covered nodes
        C = walk(GT, u, seen)  # Don't go "backward" (seen)
        seen.update(C)         # We've now seen C
        sccs.append(C)         # Another SCC found
    return sccs


TT = {
    'A': ['B', 'C'],
    'B': ['D','E'],
    'C': [],
    'D': [],
    'E': []
}
def bfs_print_tree(tree,r):
    """
    A modified single list solution
    """
    Q = []
    idx = 1
    print(str(idx)+': '+r)
    Q.extend(tree[r])
    while Q:
        idx += 1
        print(str(idx)+':')
        for u in range(0,len(Q)):
            e = Q.pop(0)
            print(e),
            Q.extend(tree[e])
        print('')

def bfs_dict(tree,r):
    """
    Bob's suggested dictionary based solution
    """
    D = {}
    idx = 1
    D[idx] = [r]
    while D[idx]:
        idx += 1
        D[idx] = []
        for u in D[idx-1]: D[idx].extend(tree[u])
    D.pop(idx) # the last dictionary element is empty - must go
    for idx in D: print(str(idx)+': '+' '.join(D[idx]))


def cloneBfree(G):
    D = {}
    for v in G:
        D[v] = {}
        for u in G[v]:            
            if not (len(G[v][u].intersection(set([(0,1)]))) == 0):
                D[v][u] = set([(0,1)])
    return D

def clrbi(G):    
    for v in G:
        d = []
        for u in G[v]:
            try:
                G[v][u].remove((edge_type['bidirected'],0))
                if len(G[v][u]) == 0:
                    d.append(u)
            except KeyError:
                pass
        for e in d:
            G[v].pop(e)

def ecj(G: dict, s: str, sccs: set = set()) -> set:
    """
    Find elementary circuits in a graph using Johnson's algorithm.
    """
    blocked = {v:False for v in G} # unblock all
    B = {v:[] for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]: unblock(w)

    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                #print stack
                sccs.add(len(stack))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s,stack)
    return sccs

def ecj_loops(G,s,sl=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G} # unblock all
    B = {v:[] for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            B[u].remove(w)
            if blocked[w]: unblock(w)
    def circuit(v, stack):
        f = False
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                f = True
                #print np.sort(stack)
                sl.add(tuple(np.sort(stack)))
            elif not blocked[u]:
                if circuit(u, stack):
                    f = True
        if f:
            unblock(v)
        else:
            for w in G[v]:
                if v not in B[w]:
                    B[w].append(v)
        stack.pop()
        return f
    circuit(s,stack)
    return sl

def gcd(a, b):
    while b != 0: a, b = b, a%b
    return a

def listgcd(l):
    if len(l)>0:
        return gcd(l[0],listgcd(l[1:]))
    else:
        return 0
def lcm(a, b): return a*b/gcd(a,b)
def chmatch(n,m,delta):
    m,n = np.sort([n,m])
    sq = np.mod(range(n,lcm(n,m)+1,n),m)
    return np.mod(delta,m) in sq

def reachable(s, G, g):
    S, Q = set(), []
    Q.append(s)
    while Q:
        u = Q.pop()
        if u in S: continue
        if g in G[u]: return True
        S.add(u)
        Q.extend(G[u])
    return False

def allpaths(G, s, g, S=[]):
    if S is None: S = []
    S.append(s)
    if s == g:
        print(S)
    else:
        for u in G[s]:
            if u in S: continue
            allpaths(G,u,g,S)
    S.remove(s)

def lz_ecj(G,s,sccs=set()):           # elementary circuit by Johnson
    blocked = {v:False for v in G}    # unblock all
    B = {v:set() for v in G}
    stack = []
    def unblock(u):
        blocked[u] = False
        for w in B[u]:
            if blocked[w]: unblock(w)
        B[u].clear()
    def circuit(v, stack):
        stack.append(v)
        blocked[v] = True
        for u in G[v]:
            if u == s:
                print('bottom')
                unblock(v)
                yield len(stack)
            elif not blocked[u]:
                print('recurse')
                for x in circuit(u, stack):
                    unblock(v)
                    yield x
            else:
                print('unmet')
                for w in G[v]: B[w].add(v)
        stack.pop()
#    circuit(s,stack)
    for v in circuit(s,stack): yield v

def iterate_allpaths(G, s, g, d=0, S=[], c=True):
    if S is None: S = []
    S.append(s)
    d += 1
    if s == g:
        if c:
            yield d-1
        else:
            yield list(S)
    else:
        for u in G[s]:
            if u in S: continue
            for v in iterate_allpaths(G,u,g,d,S,c):
                yield v
    S.remove(s)

def iddfs(G, s, g): # iterative depening DFS paths
    yielded = set()
    def recurse(G, s, g, d, S=None):
        if s not in yielded:
            yielded.add(s)
        if d == 0: return
        if S is None: S = []
        S.append(s)
        if s == g:
            yield list(S)
        else:
            for u in G[s]:
                if u in S: continue
                for v in recurse(G, u, g, d-1, S):
                    yield v
        S.remove(s)
    n = len(G)
    for d in range(n):
        #if len(yielded) == n: break
        for u in recurse(G, s, g, d):
            yield u

def reached_at_step(G, s, d):
    """
    Iterate over nodes in G reachable from s in exactly d steps
    """
    yielded = set()
    def recurse(G, s, d, B=None):
        if d == 0:
            if s not in yielded: # this avoids yielding duplicates
                yielded.add(s)
                yield s
            return
        if B is None: B = [] # black - backed out of this path
        for u in G[s]:
            #if u in B: continue
            if G[s][u] == (edge_type['bidirected'],0): continue
            for v in recurse(G, u, d-1, B):
                yield v
        B.append(s)
    for u in recurse(G, s, d):
            yield u

def d_trek(h,G,a,b,d):
    """
    Does there exist a trek with head h connecting a and b in d steps.
    """
    return set([a,b]).issubset(reached_at_step(G,h,d))

def d_biegde(G,a,b,d):
    """
    Do  a and  b  become connected  by  a bidirectional  edge after  d
    undersamples
    """
    for i in range(1,d+1):
        for u in G:
            if d_trek(u,G,a,b,i):
                return True
    return False

def undersample(G: dict, d: int, bid: bool = True) -> dict:
    """
    Perform undersampling on a graph.
    """
    N = {}
    for u in G: 
        N.update({u:{v:set([(0,1)]) for v in reached_at_step(G,u,d+1)}})
    if bid:
        items = G.keys()
        for i in range(len(items)):
            for j in range(i+1, len(items)):
                u,v = items[i],items[j]
                if d_biegde(G,u,v,d):
                    try:
                        N[u][v].add((edge_type['bidirected'],0))
                    except KeyError:
                        N[u].update({v:set([(edge_type['bidirected'],0)])})
                    try:
                        N[v][u].add((edge_type['bidirected'],0))
                    except KeyError:
                        N[v].update({u:set([(edge_type['bidirected'],0)])})
    return N

import itertools as itt
def exist_equal_paths(h,G,a,b):
    Sa, Sb, Dff = set(), set(), set()
    ag=iterate_allpaths(G,h,a,0,[],True)
    bg=iterate_allpaths(G,h,b,0,[],True)
    for v in itt.zip_longest(ag,bg):
        print(v)
        Sa.add(v[0])
        Sb.add(v[1])
        if v[0] in Sb or v[1] in Sa: return True
    return False

# checks if there  exist exact length paths from the  head node to the
# nodes at question  by iterative deepining to avoid  oing through all
# paths
def iexist_equal_paths(h,G,a,b):
    Sa, Sb = set(), set()
    Pa, Pb = [], []
    ag=iddfs(G,h,a)
    bg=iddfs(G,h,b)
    for v in itt.izip(ag,bg):
        print(v)
        Sa.add(len(v[0])); Pa.append(v[0])
        Sb.add(len(v[1])); Pb.append(v[1])
        if len(v[0]) in Sb or len(v[1]) in Sa: return True
    return False

# check  if two  unequal  length  paths can  be  compensated by  their
# elementary cycles

def has_unit_cycle(G,path):
    for v in path:
        if v in G[v]: return True
    return False
def ecj_compat(G,p1,p2):
    n = len(p1)
    m = len(p2)
    p2, p1 = [[p1,p2][i] for i in np.argsort([n,m])]
    m,n = np.sort([n,m])
    delta = n - m
    if not delta: return True # redundant check for 0
    if has_unit_cycle(G,p2): return True # equivalent
    # if the shorter path does not have cycles they are not compatible
    # if the  longer path does not  have cycles: check  if the shorter
    #                                            path    has    cycles
    #                                            divisible by delta

    # otherwise start checking
    print(p1, p2, n, m, delta)


def wc(n):
    n = n*3
    a={str(v):set([str(v+1),str(v+2),str(v+3)]) for v in range(1,n,3)}
    b={str(v):set([str(v+2)]) for v in range(2,n,3)}
    c={str(v):set([str(v+1)]) for v in range(3,n,3)}
    a.update(b)
    a.update(c)
    a.update({str(n):set()})
    return a

# Frobenius number from here: http://cgi.gladman.plus.com/wp/?page_id=563 
def residue_table(a):
  n = [0] + [None] * (a[0] - 1)
  for i in range(1, len(a)):
    d = gcd(a[0], a[i])
    for r in range(d):
      try:
        nn = min(n[q] for q in range(r, a[0], d) if n[q] != None)
      except:
        continue
      if nn != None: 
        for c in range(a[0] // d):
          nn += a[i]
          p = nn % a[0]
          nn = min(nn, n[p]) if n[p] != None else nn
          n[p] = nn
  return n
 
def frobenius_number(a):
  return max(residue_table(sorted(a))) - min(a)
 
def isSclique(G):
    n = len(G)
    for v in G:
        if sum([(0,1) in G[v][w] for w in G[v]]) < n: return False
        if sum([(2,0) in G[v][w] for w in G[v]]) < n-1: return False
    return True

# Jianyu does not use bidirected edges
def isJclique(G):
    return (sum([len(G[w].keys()) for w in G]) == len(G)**2)

def directed_inc(G,D):
    G_un = {}
    # directed edges
    for v in D:
        G_un[v] = {}
        for w in [el for el in D[v] if (0,1) in D[v][el]]:
            for e in G[w]:
                G_un[v][e] = set([(0,1)])
    return G_un
def bidirected_inc(G,D):
    # bidirected edges
    for w in D:
        # transfer old bidirected edges
        l = [e for e in D[w] if (2,0) in D[w][e]]
        for p in l:
            try: 
                 G[w][p].add((2,0))
            except KeyError: 
                 G[w][p] = set([(2,0)])
        # new bidirected dges
        l = [e for e in D[w] if (0,1) in D[w][e]]
        for p in list(combinations(l, 2)):
            try: 
                 G[p[0]][p[1]].add((2,0))
            except KeyError: 
                 G[p[0]][p[1]] = set([(2,0)])
            try: 
                G[p[1]][p[0]].add((2,0))
            except KeyError: 
                G[p[1]][p[0]] = set([(2,0)])
    return G
def increment_u(G_star, G_u):
    # directed edges
    G_un = directed_inc(G_star,G_u)
    # bidirected edges
    G_un = bidirected_inc(G_un,G_u)
    return G_un

def sample_graph(graph_g,steps=5):
    graph_g_list = [graph_g]
    for i in range(0,steps):
        g = increment_u(graph_g,graph_g_list[-1])
        graph_g_list.append(g)
    return graph_g_list
