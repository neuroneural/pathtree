import networkx
from bfutils import g2num, num2CG
from comparison import graph2nx
from ecj import undersample
from numpy import argsort

def simple_loops(g, u):
    """
    iterator over the list of simple loops of graph g at the undersample rate u
    """
    gx = graph2nx(num2CG(g2num(undersample(g,u)), len(g)))
    for l in networkx.simple_cycles(gx):
        yield l

def print_loops(g, u):
    l = [x for x in simple_loops(g,u)]
    lens = map(len, l)    
    for i in argsort(lens):
        print l[i]

def ul(l):
    """
    returns list elements that are present only once
    """
    u, r = set(), set()
    for e in l:
        if e not in u:
            u.add(e)
        else:
            r.add(e)
    return u.difference(r)
