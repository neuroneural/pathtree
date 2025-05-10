import sys

sys.path.append('./tools/')

import bclique as bq
import pathtreetools as ptt
from pathtree import PathTree
from pathtreetools import forest_to_set, to_path_forest, refine_edge
from ortools.constraint_solver import pywrapcp


class SolutionNotFoundInTime(Exception):
    """Raised when the algorithm cannot extend the PathTree within a given cutoff."""
    pass
def can_add_loop(pt, num, elements):
    """
    Check if a loop can be added to a PathTree while maintaining compatibility.
    :param pt: PathTree object
    :param num: Integer representing the loop to be added
    :param elements: List of elements for compatibility checking
    :return: Boolean indicating if the loop can be added
    """
    r = False

    # Check if the loop itself is compatible
    if ptt.isptelement(PathTree({num - pt.preset}, pre=pt.preset), num):
        r = True

    # Recursively check loopsets in the PathTree
    for e in pt.loopset:
        if isinstance(e, int):
            if can_add_loop(PathTree({e}, pre=pt.preset), num, elements):
                r = True
        elif isinstance(e, PathTree):
            if can_add_loop(e, num, elements):
                r = True

    return r


def learn_path_tree(pt, target_lags):
    """
    Return a PathTree that represents the union of observed delays.
    If target_lags is already a PathTree (i.e. pt already has the refined structure, including children),
    simply return it; otherwise, build a new tree using the observed delays.
    """
    # If target_lags is a PathTree, assume pt is already refined.
    if isinstance(target_lags, PathTree):
        return pt  # Preserve children, self-loop expansion, etc.
    
    # Otherwise, assume target_lags is a set of observed delays.
    observed = target_lags
    if not observed:
        raise ValueError("Target lag set is empty.")
    
    # Optionally, sort and check representability (if needed)
    elements = sorted(list(observed), key=lambda x: x if isinstance(x, int) else float('inf'))
    newpt = PathTree(preset=observed)
    not_representable = []
    for candidate in elements:
        if not ptt.isptelement(newpt, candidate, maxloop=10 * len(elements)):
            not_representable.append(candidate)
    if not_representable:
        print(f"Warning: the following candidates are not representable in the PathTree: {not_representable}")
    
    return newpt


def convert_bclique(bc):
    # Assume bc is a set of tuples (u, v).
    V1 = set()
    V2 = set()
    for (u, v) in bc:
        V1.add(u)
        V2.add(v)
    # Return frozensets instead of regular sets.
    return (frozenset(V1), frozenset(V2))

def bpts(graph, bclique, max_new_loops=2, maxloop=10):
    """
    Batch-PTS: refine all edges in one b-clique simultaneously,
    sharing up to `max_new_loops` new latent loops L[k].
    Returns a dict mapping (u,v) → refined PathTree.
    """
    solver = pywrapcp.Solver("batch_pts")
    solver.set_time_limit(30_000)
    V1, V2 = bclique
    edge_lag_dict: dict[tuple,str] = {}
    for u in V1:
        for v in V2:
            raw = graph[u][v].get(1) or graph[u][v].get(2)
            forest = to_path_forest(raw)
            edge_lag_dict[(u,v)] = forest_to_set(forest)
    edges = list(edge_lag_dict.items())  # [ ((u,v), E_uv), ... ]

    # 1) Shared new‐loop vars
    L = [solver.IntVar(1, maxloop, f"L[{k}]") for k in range(max_new_loops)]
    #L = [solver.IntVar(1, maxloop) for k in range(max_new_loops)]
    # 2) Per‐edge boolean weights for existing & new loops
    w_ex, w_new = {}, {}
    for (u,v), E in edges:
        forest0 = to_path_forest(graph[u][v].get(1) or graph[u][v].get(2))
        existing = [ (next(iter(pt.preset)) if isinstance(pt.preset, set)
                       else pt.preset)
                     for pt in forest0 ]

        w_ex[(u,v)]  = [solver.IntVar(0,1, f"w_ex_{u}_{v}_{i}") 
                         for i in range(len(existing))]
        # w_ex[(u,v)]  = [solver.IntVar(0,1) 
        #                  for i in range(len(existing))]
        w_new[(u,v)] = [solver.IntVar(0,1, f"w_new_{u}_{v}_{k}") 
                         for k in range(max_new_loops)]
        # w_new[(u,v)] = [solver.IntVar(0,1) 
        #                  for k in range(max_new_loops)]

        # 3) build linearized “product” for new loops: p_k = L[k] * w_new[k]
        M = maxloop
        p_vars = []
        for k, wv in enumerate(w_new[(u,v)]):
            pk = solver.IntVar(0, M, f"p_{u}_{v}_{k}")
            #pk = solver.IntVar(0, M)
            solver.Add(pk <= L[k])
            solver.Add(pk <= wv * M)
            solver.Add(pk >= L[k] - (1 - wv) * M)
            p_vars.append(pk)

        # 4) sum up existing & new‐loop contributions
        all_vars = w_ex[(u,v)] + p_vars
        coeffs   = existing + [1]*len(p_vars)
        sum_var  = solver.ScalProd(all_vars, coeffs)

        # 5) enforce (sum_var + preset) ∈ E exactly
        solver.Add(solver.MemberCt(sum_var + min(E), list(E)))

    # 6) solve once
    all_vars = L + sum(w_ex.values(), []) + sum(w_new.values(), [])
    solution = solver.Assignment()
    solution.Add(all_vars)
    db = solver.Phase(all_vars,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)
    if not solver.NextSolution():
        solver.EndSearch()
        raise RuntimeError("No solution for b-clique")
    # capture
    L_vals     = [solution.Value(v)       for v in L]
    w_ex_vals  = {e:[solution.Value(wv) for wv in w_ex[e]]  for e in w_ex}
    w_new_vals = {e:[solution.Value(wv) for wv in w_new[e]] for e in w_new}
    solver.EndSearch()

    # 7) rebuild & refine a PathTree for each edge
    out = {}
    for (u,v), E in edges:
        forest0 = to_path_forest(graph[u][v].get(1) or graph[u][v].get(2))
        preset  = min(E)

        # keep existing loops if solver didn’t turn them off
        existing = [ (next(iter(pt.preset)) if isinstance(pt.preset, set)
                        else pt.preset)
                     for pt in forest0 ]
        loops = { coef for coef, flag in zip(existing, w_ex_vals[(u,v)]) if flag == 0 }
        # add any new loops the solver activated
        loops |= { L_vals[k] for k, flag in enumerate(w_new_vals[(u,v)]) if flag == 1 }

        # 2‐stage: build the raw tree, then *greedily* refine to cover *all* E
        raw_pt  = PathTree(preset=preset, loopset=loops)
        full_pt = refine_edge(raw_pt, E)

        out[(u,v)] = full_pt
    return list(out.values())

def getagenerator(g):
    """
    Generate b-cliques and compute associated PathTrees.
    :param g: Input graph
    :return: List of b-cliques
    """
    bclqs = bq.bcliques(g)
    for clq in bclqs:
        bpts(clq)
    return bclqs
