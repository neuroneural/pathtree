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
    Extend the given PathTree pt so that it can represent the observed edge-lag set target_lags.
    Instead of iterating over all integers, we only consider candidate values from target_lags.
    If a candidate cannot be represented by the current tree, it is simply skipped.
    """
    # Instead of taking just the minimum delay, store the entire set.
    if not target_lags:
        raise ValueError("Target lag set is empty.")
    newpt = PathTree(preset=target_lags)

    def rpath(remaining, npt):
        for candidate in remaining:
            # If npt already generates candidate, nothing to do.
            if ptt.isptelement(npt, candidate, maxloop=10 * len(elements)):
                continue
            # Try to extend npt to generate candidate.
            try:
                refined = growtree(npt, candidate, elements)
            except SolutionNotFoundInTime:
                # Skip candidate if no extension is found.
                continue
            if refined is None:
                # Candidate not representable; skip it.
                continue
            # Otherwise, update the current tree.
            npt = refined
        return npt

    # Process candidates (other than the smallest, which is already the preset)
    refined_pt = rpath(elements[1:], newpt)

    # Instead of throwing an exception if some candidate is not representable,
    # we log a warning and return the tree as-is.
    representable = []
    not_representable = []
    for candidate in elements:
        if ptt.isptelement(refined_pt, candidate, maxloop=10 * len(elements)):
            representable.append(candidate)
        else:
            not_representable.append(candidate)
    if not_representable:
        print(f"Warning: the following elements are not representable in the PathTree: {not_representable}")
    return refined_pt

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
