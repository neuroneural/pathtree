import sys

sys.path.append('./tools/')
import graphkit as gk
import numpy as np
from ortools.constraint_solver import pywrapcp

def grelabel(g):
    """
    Relabel graph nodes to numerical contiguous labels 1 through len(g)
    :param g: input graph
    :return: relabeled graph and name mapping
    """
    gg = {}
    nodemap = {x[0]: x[1] + 1 for x in zip(list(g.keys()), range(len(g)))}
    for v in g:
        gg[nodemap[v]] = {}
        for w in g[v]:
            gg[nodemap[v]][nodemap[w]] = g[v][w]
    return gg, {v: k for k, v in nodemap.items()}


def crelabel(c, map):
    newcliq = set()
    for e in c:
        newcliq.add((map[e[0] + 1], map[e[1] + 1]))
    return newcliq


def printedges(g):
    l = gk.edgelist(g)
    for e in l:
        print(e[0] - 1, '->', e[1] - 1)


def edge_exists(i, j, g):
    return j + 1 in g[i + 1]


def numkids(i, edges, solver):
    N = len(edges)
    return solver.Sum([edges[i][k] for k in range(N)])


def numparents(i, edges, solver):
    N = len(edges)
    return solver.Sum([edges[k][i] for k in range(N)])


def clique_constrain(solver, parents, children, edges, g):
    N = len(parents)

    # declare constraints
    solver.Add(solver.Sum(parents) > 0)
    solver.Add(solver.Sum(children) > 0)

    for i in range(N):
        # makes sure that there exists at least one outgoing edge from node i if it is marked in parents for this b-clique
        solver.Add((parents[i] == 1) == (solver.Sum([edges[i][k] for k in range(N)]) >= 1))
        # this makes sure that there exists at least one incoming edge to node i if it is marked as a child
        solver.Add((children[i] == 1) == (solver.Sum([edges[k][i] for k in range(N)]) >= 1))

        solver.Add((children[i] == 1) == (solver.Sum(parents) <= numparents(i, edges, solver)))
        solver.Add((parents[i] == 1) == (solver.Sum(children) <= numkids(i, edges, solver)))

        for j in range(N):
            # edge existence constraints
            if not edge_exists(i, j, g):
                solver.Add(edges[i][j] == 0)
            else:
                solver.Add((edges[i][j] == 1) == (parents[i] * children[j] == 1))


def bcliques(g, verbose=False):
    solver = pywrapcp.Solver("b-clique")
    g, mp = grelabel(g)
    # declare variables
    edges = []
    N = len(g)
    for i in range(N):
        e = []
        for j in range(N):
            e.append(solver.IntVar(0, 1, "%i -> %i" % (i, j)))
        edges.append(e)

    parents = [solver.IntVar(0, 1, "%i" % (i)) for i in range(N)]
    children = [solver.IntVar(0, 1, "%i" % (i)) for i in range(N)]

    # declare constraints
    clique_constrain(solver, parents, children, edges, g)

    # run the solver
    solution = solver.Assignment()
    solution.Add([edges[i][j] for i in range(N) for j in range(N)])
    solution.Add(parents)
    solution.Add(children)

    collector = solver.AllSolutionCollector(solution)
    solver.Solve(solver.Phase([edges[i][j] for i in range(N) for j in range(N)] + children + parents,
                              solver.CHOOSE_FIRST_UNBOUND,
                              solver.ASSIGN_MAX_VALUE),
                 [collector])
    num_solutions = collector.SolutionCount()

    # output solutions
    if verbose:
        print("num_solutions:", num_solutions)
        print("failures:", solver.Failures())
        print("branches:", solver.Branches())
        print("WallTime:", solver.WallTime())

    cliques = []
    pts = set()
    for s in range(num_solutions):
        c = set()
        qval = [collector.Value(s, edges[i][j]) for i in range(N) for j in range(N)]
        pval = [collector.Value(s, parents[i]) for i in range(N)]
        cval = [collector.Value(s, children[i]) for i in range(N)]
        if verbose:
            print(" ----------------- ", s)
        if verbose:
            print(np.asarray([pval, cval]))
        for i in range(len(qval)):
            if qval[i]:
                e = np.unravel_index(i, [N, N])
                c.add((e[0], e[1]))
                if verbose:
                    print(e[0], "->", e[1])
        if not (True in list(map(lambda x: c.issubset(x), cliques))):
            if not c.issubset(pts):
                pts = pts.union(c)
                cliques.append(c)
        if verbose:
            print("---")

    # check for b-cliques for which all edges are covered in other b-cliques
    cc = []
    for i in range(len(cliques)):
        bcl = cliques.pop()
        pts = set()
        for j in range(len(cliques)):
            pts = pts.union(cliques[j])
        for j in range(len(cc)):
            pts = pts.union(cc[j])

        if not bcl.issubset(pts):
            cc.append(bcl)
    return list(map(lambda x: crelabel(x, mp), cc))
