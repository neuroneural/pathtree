import sys

sys.path.append('./tools/')

import bclique as bq
import latents as lt
import pathtreetools as ptt
from pathtree import PathTree
from pathtreetools import compute_edge_lags, growtree


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

def bpts(graph, bclique):
    # Convert bc to a tuple (V1, V2)
    bc_converted = convert_bclique(bclique)
    edge_lag_dict = compute_edge_lags(graph, bc_converted)
    pts = set()
    for edge in edge_lag_dict:
        target_lags = edge_lag_dict[edge]
        # If the target lag is already a PathTree, use it; otherwise, create one.
        if isinstance(target_lags, PathTree):
            pt = target_lags
        else:
            pt = PathTree(preset=min(target_lags, key=lambda x: x if isinstance(x, int) else float('inf')))
        try:
            learned_pt = learn_path_tree(pt, target_lags)
            pts.add(learned_pt)
        except ValueError as e:
            print(e)
            continue
    return pts

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
