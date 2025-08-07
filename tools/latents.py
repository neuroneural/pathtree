import sys, os

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np
from random import shuffle

from copy import deepcopy
from pathtree import PathTree, osumset_full
from pathtreetools import to_path_forest, forest_to_set

g_var_counter = 0

def get_fresh_loop_var():
    """
    Returns a fresh variable name each call, i.e. 'k', 'l', 'm', 'n', ...
    Once it reaches 'z', you may want a more robust approach, but this is 
    enough for demonstration.
    """
    global g_var_counter
    # 'a' is chr(97). Next call is chr(98)='b', then 'c', ...
    var_name = chr(97 + g_var_counter)
    g_var_counter += 1
    return var_name

# ---------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------
def _is_helper(node: int, g: dict) -> bool:
    """
    Return True  ⇢  node is a 1-lag relay that decompress_to_unit_graph
                    created (one 1-lag tag-1 in-edge, one 1-lag tag-1
                    out-edge, no self-loop, no tag-2).
           False ⇢  anything else (observed vertices, bidirected latents,
                    A-P compression latents, etc.).
    """
    if node not in g:
        return False

    # ----- outgoing -----
    out = g[node]
    if len(out) != 1:
        return False
    (child, edict_out), = out.items()
    if edict_out != {1: {1}}:
        return False

    # ----- incoming -----
    incoming = [
        edict_in
        for par, nbrs in g.items()
        if node in nbrs
        for edict_in in [nbrs[node]]
    ]
    if len(incoming) != 1 or incoming[0] != {1: {1}}:
        return False

    # ----- self-loop or bidirected? -----
    if node in out and 2 in out[node]:
        return False

    return True

def update_directed_edges(graph, r):
    """
    Update directed edges when removing a latent node 'r'.
    This handles recalculating lags between parent and child nodes of 'r',
    correctly accounting for self-loops.

    :param graph: The graph as a dictionary.
    :param r: The node to be removed.
    :return: The updated graph.
    """
    parents_r = {p: graph[p][r][1] for p in graph if r in graph[p]}
    children_r = {c: graph[r][c] for c in graph[r]}

    # Detect self-loop at r
    self_loop = graph[r][r][1] if r in graph[r] else set()
    print(f"Self-loop at {r}: {self_loop}")

    # Expand self-loop into infinite set approximation
    if self_loop:
        self_loop_expanded = set()
        for x in self_loop:
            self_loop_expanded.update({x * k for k in range(10)})  # Approximate infinite sequence
    else:
        self_loop_expanded = {0}  # No self-loop, just pass 0
    print(f"Expanded self-loop contributions: {self_loop_expanded}")

    updated_graph = deepcopy(graph)
    del updated_graph[r]  # Remove the latent node

    for p in parents_r:
        for c in children_r:
            # Combine parent paths, self-loop repetitions, and child paths
            combined_lags = osumset_full(parents_r[p], self_loop_expanded)
            combined_lags = osumset_full(combined_lags, children_r[c].get(1, set()))
            print(f"Updated edge {p} -> {c} lags: {combined_lags}")

            if p == c:  # Self-loop case
                if p in updated_graph and p in updated_graph[p]:
                    updated_graph[p][p][1].update(combined_lags)
                else:
                    updated_graph.setdefault(p, {})[p] = {1: combined_lags}
            else:  # Regular edge case
                if c in updated_graph[p]:
                    updated_graph[p][c][1].update(combined_lags)
                else:
                    updated_graph[p][c] = {1: combined_lags}

    # Clean up remaining references to r in the graph
    for node in list(updated_graph.keys()):
        updated_graph[node].pop(r, None)

    print("Updated graph after removing node:", updated_graph)
    return updated_graph

def bidirected_forest(ah_forest, hb_forest):
    """
    For induced bidirected edges, compute (y - x) for every base in ah_forest and hb_forest,
    returning a forest (list of PathTrees) with each difference as the preset.
    """
    new_forest = []
    for pt1 in ah_forest:
        bases1 = get_base_as_set(pt1.preset)
        for pt2 in hb_forest:
            bases2 = get_base_as_set(pt2.preset)
            for x in bases1:
                for y in bases2:
                    diff_val = abs(y - x)
                    new_forest.append(PathTree(preset={diff_val}, loopset=set()))
    return new_forest

def update_bidirected_edges(graph, r):
    """
    Update bi-directed edges when removing a latent node 'r'.
    Recalculate lags between parents and children of 'r' by taking the
    difference in lags, then store the final edge under whichever node
    has the smaller lag (the "closer" node), pointing to the "further" node.

    This ensures that if node 1->2 has lag=6 and node 1->3 has lag=7,
    node 2 is closer (6 < 7), so we store final '2 <-> 3 = (1)', not '3 <-> 2'.
    """

    # If 'r' doesn't exist, just return
    if r not in graph:
        print(f"Skipping Node {r}: It does not exist in the graph.")
        return graph

    # Identify which nodes connect bidirectionally to 'r'
    parents_r = {
        p: graph[p][r].get(2, set())
        for p in graph
        if (r in graph[p]) and (2 in graph[p][r])
    }
    children_r = {
        c: graph[r][c].get(2, set())
        for c in graph.get(r, {})
    }

    print(f"Node to Remove: {r}")
    print(f"Parents of {r} (Bi-directed): {parents_r}")
    print(f"Children of {r} (Bi-directed): {children_r}")

    updated_graph = deepcopy(graph)

    # Remove 'r' itself from the adjacency
    if r in updated_graph:
        del updated_graph[r]

    # For each (p, c) pair, we produce a correlation p <-> c
    # by taking the difference in lags: abs(lp - lc)
    for p in parents_r:
        for c in children_r:
            # Collect all differences
            diff_lags = {
                abs(lp - lc)
                for lp in parents_r[p]
                for lc in children_r[c]
            }

            if p == c:
                # Self-loop scenario
                if p in updated_graph and p in updated_graph[p]:
                    updated_graph[p][p].setdefault(2, set()).update(diff_lags)
                else:
                    updated_graph.setdefault(p, {})[p] = {2: diff_lags}
                continue

            # Decide which node is "closer" to r
            # We'll compare the MAX of each node's lag set. Whichever is smaller is 'closer'.
            max_p = max(parents_r[p]) if parents_r[p] else 0
            max_c = max(children_r[c]) if children_r[c] else 0

            if max_p < max_c:
                # p is closer => store p->c
                closer, further = p, c
            elif max_c < max_p:
                # c is closer or they tie => store c->p
                closer, further = c, p
            else:
                # tie — pick the lexicographically smaller node as “closer”
                if str(p) < str(c):
                    closer, further = p, c
                else:
                    closer, further = c, p
            # Insert the final bidirected edge in updated_graph[closer][further][2]
            if closer not in updated_graph:
                updated_graph[closer] = {}
            if further not in updated_graph[closer]:
                updated_graph[closer][further] = {}
            updated_graph[closer][further].setdefault(2, set()).update(diff_lags)

    # Clean up references to r
    for node in list(updated_graph.keys()):
        if r in updated_graph[node]:
            del updated_graph[node][r]

    return updated_graph

def sequential_forward_inference(graph, latent_nodes):
    """
    Perform sequential forward inference by marginalizing over latent nodes.

    :param graph: The input graph as a dictionary.
    :param latent_nodes: A set of latent nodes to marginalize.
    :return: The updated graph.
    """
    updated_graph = deepcopy(graph)

    for r in latent_nodes:
        print(f"Removing Node {r}...")
        if r not in updated_graph:
            print(f"Skipping Node {r}: Already removed.")
            continue

        # Update directed and bi-directed edges
        updated_graph = update_directed_edges(updated_graph, r)
        updated_graph = update_bidirected_edges(updated_graph, r)

        # Remove the node `r` itself
        if r in updated_graph:
            del updated_graph[r]

        # Clean up any references to `r`
        for node in list(updated_graph.keys()):
            if r in updated_graph[node]:
                del updated_graph[node][r]

    return updated_graph


def g2lg(g):
    """
    Convert a data structure encoding the MSL-type graph into a structure encoding latents graph
    :return: a graph with integer indices and sets as weights

    Args:
        g (MSL-graph): 
    """
    edge_type = {(0, 1): 1, (2, 0): 2}
    edge_weight = {(0, 1): 1, (2, 0): 0}
    lg = {int(e): {int(c): {edge_type[w]: {edge_weight[w]}
                            for w in g[e][c]}
                   for c in g[e]}
          for e in g}
    return fix_selfloops(lg)


def fix_selfloops(g):
    for v in g:
        if v in g[v]:
            g[v][v] = {1: {PathTree({1})}}
    return g


def gtranspose(G):  # Transpose (rev. edges of) G
    GT = {u: {} for u in G}
    for u in G:
        for v in G[u]:
            if 1 in G[u][v]:
                GT[v][u] = {1: {1}}  # Add all reverse edges
    return GT


def parents(g, N):  # an inefficient way to do it
    t = gtranspose(g)  # Transpose the graph
    return {n: g[n][N][1]
            for n in t[N] if n != N}  # Return all children

def children(g, N):
    """Return children of N, ensuring all edge types are considered."""
    return {
        n: g[N][n][1] if 1 in g[N][n] else g[N][n][2]
        for n in g[N] if 1 in g[N][n] or 2 in g[N][n]
    }

def remove_node(g, N):
    del g[N]
    for V in g:
        if N in g[V]:
            del g[V][N]


def iterate_ws(ws):
    starts = []
    for e in ws:
        if type(e) is int:
            starts.append(e)
            continue
        if type(e) is set:
            starts.extend(list(e))
            continue
        starts.append(e.preset)
        starts.append([e.preset, iterate_ws(e.loopset)])
    return starts


def iterate_pt(pt):  # iterate over a path tree
    starts = [pt.preset]
    for e in pt.loopset:
        if type(e) is int:
            starts.append(e)
            continue
        if type(e) is set:
            starts.extend(list(e))
            continue
        starts.append(e.preset)
        starts.append([e.preset, iterate_ws(e.loopset)])
    return starts

def get_base_as_set(x):
    """
    Recursively converts x (which may be an int, a set, a list, or a PathTree)
    into a plain set of integers.
    If x is an int, return {x}.
    If x is a list or set, flatten each element.
    If x is a PathTree, recurse on its preset.
    """
    if isinstance(x, int):
        return {x}
    elif isinstance(x, (list, set)):
        result = set()
        for item in x:
            result |= get_base_as_set(item)
        return result
    elif isinstance(x, PathTree):
        # We assume that x.preset is either an int or a set or list.
        return get_base_as_set(x.preset)
    else:
        raise TypeError(f"get_base_as_set expects an int, list, set, or PathTree, got {type(x)}")

    
def sum_forests(f1, f2):
    """
    For each PathTree in f1 and each PathTree in f2, produce new PathTrees
    whose preset is the Cartesian sum of their bases and whose loopset is the union.
    Return the new forest (list of PathTrees).
    """
    new_forest = []
    for pt1 in f1:
        bases1 = get_base_as_set(pt1.preset)
        loops1 = set(pt1.loopset)
        for pt2 in f2:
            bases2 = get_base_as_set(pt2.preset)
            loops2 = set(pt2.loopset)
            s = osumset_full(bases1, bases2)
            if 0 in s:
                s.discard(0)
            merged_loops = loops1.union(loops2)
            for b in s:
                new_forest.append(PathTree(preset={b}, loopset=merged_loops.copy()))
    return new_forest

def unify_forests(f1, f2):
    """
    Merge each PathTree in f2 into f1 if they have the same preset (i.e. same base).
    Otherwise, keep them separate.
    """
    new_forest = f1[:]
    for pt2 in f2:
        base2 = get_base_as_set(pt2.preset)
        inserted = False
        for pt1 in new_forest:
            base1 = get_base_as_set(pt1.preset)
            if base1 == base2:
                pt1.loopset.update(pt2.loopset)
                inserted = True
                break
        if not inserted:
            new_forest.append(pt2)
    return new_forest

def merge_weightsets(ab, ah, hb, hh, induced_bidirected=False):
    """
    Combine delay contributions from parent->H (ah) and H->child (hb) and incorporate 
    self-loop delays at H (hh), returning a path-forest (a list of PathTrees) with:
      - preset = sum of delays (possibly from a direct edge if available), and 
      - loopset = union of extra delays (from indirect paths and self-loops).
    """
    # Convert ah and hb into forests.
    ah_forest = to_path_forest(ah)
    hb_forest = to_path_forest(hb)
    
    if induced_bidirected:
        new_base = {abs(y - x) for pt1 in ah_forest for x in get_base_as_set(pt1.preset)
                           for pt2 in hb_forest for y in get_base_as_set(pt2.preset)}
        if 0 in new_base and len(new_base) > 1:
            new_base.discard(0)
        pt = PathTree(preset=new_base, loopset=set())
        return [pt]
    else:
        # Instead of summing just the presets, use sum_forests so that the loopset info is preserved.
        indirect_forest = sum_forests(ah_forest, hb_forest)
        # If multiple branches result, merge them into one using unify_forests.
        # Cap the "A-P loop" optimisation: compress only when we have ≥ 3 bases
        if len(indirect_forest) == 1:
            bases = get_base_as_set(indirect_forest[0].preset)  # {int,…}
            if len(bases) <= 2:
                # Keep it flat: no <d> loop
                pt_indirect = PathTree(preset=bases, loopset=set())
            else:
                pt_indirect = indirect_forest[0]                # OK to keep loop
        else:
            # Multiple branches → unify as before (may produce a loop root)
            indirect_forest = unify_forests(indirect_forest, [])
            pt_indirect = indirect_forest[0]
        bases = get_base_as_set(pt_indirect.preset)     # {int,…}
        if len(bases) <= 2:
            pt_indirect = PathTree(preset=bases, loopset=set())
        new_preset = pt_indirect.preset
        
        # If a direct edge delay exists, add it to every branch.
        direct = get_base_as_set(ab)
        if direct:
            new_presets = set()
            for d in direct:
                for b in get_base_as_set(new_preset):
                    new_presets.add(d + b)
            new_preset = new_presets
        
        # Process self-loop delays (hh)
        extra_hh = set()     # flat self-loop delays
        nested_hh = set()    # nested PathTree delays
        if hh:
            if isinstance(hh, list):
                for item in hh:
                    if hasattr(item, 'preset'):
                        nested_hh.add(item)
                    else:
                        extra_hh = extra_hh.union(get_base_as_set(item))
            else:
                extra_hh = extra_hh.union(get_base_as_set(hh))
        
        bases = get_base_as_set(new_preset)
        if len(bases) == 1:                     # ← new guard
            extra_hh  = set()
            nested_hh = set()
        # Combine the self-loop delays with those already in the indirect forest.
        final_loopset = pt_indirect.loopset.union(extra_hh).union(nested_hh)
        pt_final = PathTree(preset=new_preset, loopset=final_loopset)
        #print("Final merged PathForest:", forest)
        return [pt_final]
def add_self_loops(forest, hh):
    """
    Attach self-loop delays from hh to each PathTree in the forest.
    Instead of overwriting the loopset, we merge them:
      - For a tree with base b, if b==2 then attach 10 if present;
      - If b==3 then attach 27 if present.
    (Adjust these heuristics as needed.)
    """
    new_forest = []
    for pt in forest:
        b = next(iter(pt.preset))  # assume preset is a singleton
        # Start with the existing loopset
        attached = set(pt.loopset)
        for L in hh:
            if b == 2 and L == 10:
                attached.add(L)
            elif b == 3 and L == 27:
                attached.add(L)
        new_forest.append(PathTree(preset=pt.preset, loopset=attached))
    return new_forest
def flatten_to_ints(obj):
    """
    Recursively flatten obj (which may be an int, a PathTree, or a set/list of them)
    into a plain set of integers. If a PathTree has preset = {6}, we add 6, etc.
    """
    results = set()

    # If obj is None or empty, return an empty set
    if not obj:
        return results

    # If obj is a single item, wrap it in a list for uniform handling
    if isinstance(obj, (int, PathTree, set, list)):
        items = obj if isinstance(obj, (set, list)) else [obj]
    else:
        # Not a known type => skip or raise error
        return results

    for item in items:
        # If int, just add it
        if isinstance(item, int):
            results.add(item)

        # If a PathTree => flatten its preset
        elif isinstance(item, PathTree):
            # item.preset might be an int or a set of ints
            if isinstance(item.preset, int):
                results.add(item.preset)
            elif isinstance(item.preset, set):
                # add them all
                for val in item.preset:
                    if isinstance(val, int):
                        results.add(val)
                    # If val is another PathTree, you could even go deeper
            # If you want to handle item.loopset, do so as well
            # but if it's empty, that's fine.

        # If item is a set => flatten that recursively
        elif isinstance(item, set):
            results |= flatten_to_ints(item)

        # If item is a list => flatten that recursively
        elif isinstance(item, list):
            results |= flatten_to_ints(item)

        # else skip or raise TypeError
    return results

def _condense_flat_roots(forest):
    """
    Merge every *flat* root (loopset == ∅) into ONE PathTree whose preset
    is the union of all their presets.  Any non-flat roots are kept unchanged.
    """
    flat_bases = set()
    non_flat   = []
    for pt in forest:
        if pt.loopset:                       # has loops → keep as-is
            non_flat.append(pt)
        else:                                # completely flat
            flat_bases |= get_base_as_set(pt.preset)

    if flat_bases:
        non_flat.insert(0, PathTree(preset=flat_bases, loopset=set()))
    return non_flat

def merge_forest_by_base(forest):
    """
    Merge PathTrees that share the same *base* delay by taking the union of
    their loop-sets, then condense any remaining flat roots so that
    [{0},{1}] becomes [{0,1}].
    """
    merged: dict[int, set[int]] = {}             # base → union(loopset)

    for pt in forest:
        bases = get_base_as_set(pt.preset)       # {int,…}
        loops = flatten_to_ints(pt.loopset)      # {int,…}

        for b in bases:                          # (we assume singleton bases)
            merged.setdefault(b, set()).update(loops)

    # rebuild as PathTree objects
    flat_forest = [
        PathTree(preset={b}, loopset=ls) for b, ls in merged.items()
    ]

    # 🔸 final visual condensation
    return _condense_flat_roots(flat_forest)

def hide_node(g, H, *, zero_incoming=False):
    """
    Remove vertex H from graph g.
        • if zero_incoming is True the parent->H delay is forced to 0
        (used for helper relays inserted by decompress_to_unit_graph).

        • otherwise the real delays into H are preserved (ordinary latents).
    1. For each parent p of H and each child c of H, update directed edges p -> c 
       using merge_weightsets.
    2. Then, regardless of whether H has observed parents, induce bidirected edges among
       every pair (c1, c2) of children of H. For each pair, compute the induced bidirected
       edge using merge_weightsets with induced_bidirected=True, then use flatten_to_ints
       to decide which node is 'closer' (i.e. has the smaller maximum delay when flattened)
       so that the edge is stored as: closer -> further under key 2.
    3. Clean up all references to H.
    """
    gg = deepcopy(g)
    # ensure any existing bidirected lag-sets become a list of PathTree
    for u in list(gg):
        for v in list(gg[u]):   
            for et in (1, 2):
                if et in gg[u][v] and isinstance(gg[u][v][et], set):
                    gg[u][v][et] = to_path_forest(gg[u][v][et])
    if H not in g:
        raise KeyError(f"Node {H} not found in graph.")

    # Children of H: dictionary {child: delay from H->child}
    ch = {c: to_path_forest(raw) for c, raw in children(g, H).items()}
    # Parents of H: include both directed (etype=1) and bidirected (etype=2) into H
    pa_orig = {}
    for p, nbrs in g.items():
        if H in nbrs:
            for et, raw in nbrs[H].items():
                # raw may be a set(int), a PathTree, or a list[…]
                pa_orig.setdefault(p, []).extend(
                    to_path_forest(raw)          # always returns list[PathTree]
                )
    pa = {p: {0} for p in pa_orig} if zero_incoming else pa_orig
    # ---------- robust self-loop fetch ----------
    # 2) self-loop on H, if any
    def _self_loop(node):
        if node not in g or node not in g[node]:
            return set()
        raw = g[node][node]
        if isinstance(raw, dict):
            raw = raw.get(1, raw.get(2, set()))
        roots = raw if isinstance(raw, list) else [raw]
        return forest_to_set(roots) if not isinstance(raw, set) else raw

    sl = _self_loop(H)

    # Remove H from gg first.
    remove_node(gg, H)

    # --- Branch 1: update p → c edges -----------------------------------
    if pa:
        for p in pa:
            if p not in gg:          # parent was deleted earlier
                continue
            for c in ch:
                if c == H or p == H:
                    continue
                if p == c:
                    was_bidir = False      # a loop is never bidirected
                    key       = 1

                    existing   = gg[p][p].get(key) if p in gg[p] else None
                    old_forest = to_path_forest(existing) if existing else []

                    new_branch = merge_weightsets(
                        [], pa[p], ch[c], to_path_forest(sl) if sl else [],
                        induced_bidirected = False
                    )
                    payload = merge_forest_by_base(old_forest + new_branch)

                    gg.setdefault(p, {}).setdefault(p, {})[1] = payload
                    continue

                # is the composite edge directed or bidirected?
                was_bidir = (2 in g.get(p, {}).get(H, {}) and
                            2 in g.get(H, {}).get(c, {}))
                key = 2 if was_bidir else 1

                # what is already on p → c (if anything)
                existing    = gg[p][c].get(key) if c in gg[p] else None
                old_forest  = to_path_forest(existing) if existing is not None else []

                # weights
                pa_w = pa[p]
                ch_w = ch[c]
                sl_w = to_path_forest(sl) if sl else []

                # new branch through latent H
                new_branch = merge_weightsets([], pa_w, ch_w, sl_w,
                                            induced_bidirected=was_bidir)

                payload = merge_forest_by_base(old_forest + new_branch)

                # ── ZERO-LAG TIDYING for *both* key==1 and key==2 ──────────
                if key == 2:
                    pure_zero_before = (
                        old_forest and all(pt.preset == {0} for pt in old_forest)
                    )

                    if pure_zero_before:              # edge was exactly {0} already
                        for pt in payload:
                            pt.preset = {0}
                    else:
                        had_zero_old = any(0 in pt.preset for pt in old_forest)
                        for pt in payload:
                            if 0 in pt.preset and len(pt.preset) > 1 and not had_zero_old:
                                pt.preset.discard(0)
                # ────────────────────────────────────────────────────────────

                gg.setdefault(p, {}).setdefault(c, {})[key] = payload
    # --------------------------------------------------------------------


    # --- Branch 2: Always induce bidirected edges among all pairs of children of H ---
    if ch and (len(ch) > 1):
        children_list = list(ch.keys())
        for i in range(len(children_list)):
            for j in range(i + 1, len(children_list)):
                c1 = children_list[i]
                c2 = children_list[j]

                lag_c1 = ch[c1] if c1 in ch else set()
                lag_c2 = ch[c2] if c2 in ch else set()

                # Flatten delay sets into plain ints.
                f_c1 = flatten_to_ints(lag_c1)
                f_c2 = flatten_to_ints(lag_c2)

                # Compute the induced bidirected edge using merge_weightsets with induced_bidirected=True.
                forest = merge_weightsets(set(), lag_c1, lag_c2, sl, induced_bidirected=True)

                # Decide which child is "closer" (smaller max delay) after flattening.
                max_c1 = max(f_c1) if f_c1 else 0
                max_c2 = max(f_c2) if f_c2 else 0

                if max_c1 < max_c2:
                    closer, further = c1, c2
                elif max_c2 < max_c1:
                    closer, further = c2, c1
                else:
                    # tie — pick lexicographically smaller node as “closer”
                    if str(c1) < str(c2):
                        closer, further = c1, c2
                    else:
                        closer, further = c2, c1
                if closer not in gg:
                    gg[closer] = {}
                if further not in gg[closer]:
                    gg[closer][further] = {}

                # ── NEW: merge into whichever orientation already exists ─────────────
                if 2 in gg.get(closer, {}).get(further, {}):          # same orientation
                    old_forest   = to_path_forest(gg[closer][further][2])
                    payload      = merge_forest_by_base(old_forest + forest)

                    pure_zero_before = (
                        old_forest and all(pt.preset == {0} for pt in old_forest)
                    )
                    if pure_zero_before:
                        # keep it *exactly* {0}; discard any newcomers
                        for pt in payload:
                            pt.preset = {0}
                    else:
                        had_zero_old = any(0 in pt.preset for pt in old_forest)
                        for pt in payload:
                            if 0 in pt.preset and len(pt.preset) > 1 and not had_zero_old:
                                pt.preset.discard(0)

                    gg[closer][further][2] = payload

                elif 2 in gg.get(further, {}).get(closer, {}):               # opposite orientation
                    old_forest = to_path_forest(gg[further][closer][2])
                    payload    = merge_forest_by_base(old_forest + forest)

                    # ■ keep a pure {0} exactly as {0}
                    if all(pt.preset == {0} for pt in old_forest):
                        for pt in payload:
                            pt.preset = {0}

                    # store on the existing (opposite) edge and stop
                    for pt in payload:                       # never keep 0 with other lags
                        if 0 in pt.preset and len(pt.preset) > 1:
                            pt.preset.discard(0)
                    gg[further][closer][2] = payload
                    continue                                 # done with this child-pair

                else:                                        # edge doesn’t exist yet
                    payload = merge_forest_by_base(forest)

                # tidy up: drop a 0 that’s mixed with non-zeros
                for pt in payload:
                    if 0 in pt.preset and len(pt.preset) > 1:
                        pt.preset.discard(0)

                gg.setdefault(closer, {}).setdefault(further, {})[2] = payload
                # ────────────────────────────────────────────────────────────────────
    # Cleanup: remove any leftover references to H in the graph.
    for node in list(gg.keys()):
        if H in gg[node]:
            del gg[node][H]
    if H in gg:
        del gg[H]
    for u in list(gg):
        for v in list(gg[u]):
            for et in tuple(gg[u][v].keys()):
                raw = gg[u][v][et]
                # convert ANY raw value (set, int, tuple, PathTree, list) into List[PathTree]
                forest_list = to_path_forest(raw)
                gg[u][v][et] = forest_list
    # ─── 2) Loop‐child normalization: ensure every loop is a PathTree ─────
    def normalize_loops(pt: PathTree):
        new_loops = set()
        for child in pt.loopset:
            if isinstance(child, PathTree):
                normalize_loops(child)
                new_loops.add(child)
            else:
                wrapped = to_path_forest(child)[0]
                normalize_loops(wrapped)
                new_loops.add(wrapped)
        pt.loopset = new_loops

    for u in gg:
        for v in gg[u]:
            for et in list(gg[u][v].keys()):        # ← use .keys(), not .values()
                forest_list = gg[u][v][et]
                # make sure we have a list of PathTrees
                if not isinstance(forest_list, list):
                    forest_list = [forest_list]

                # (a) normalise child loops
                for pt in forest_list:
                    normalize_loops(pt)

                # (b) collapse flat roots that share the same base
                gg[u][v][et] = merge_forest_by_base(forest_list)

    return gg

def hide_nodes_old(g, nodelist, dosort=True):
    nodeset = set()  # make sure not to delete a node twice
    if dosort: nodelist = sortbydegree(nodelist, g)
    gg = deepcopy(g)
    for n in nodelist:
        if n in nodeset: continue
        gg = hide_node(gg, n)
        nodeset.add(n)
    return gg
def hide_nodes(g, latents, *, zero_incoming=False):
    gg        = deepcopy(g)
    remaining = set(latents)
    while remaining:
        changed = False
        for H in list(remaining):
            if H not in gg:
                remaining.discard(H)
                continue
            gg_new = hide_node(gg, H, zero_incoming=zero_incoming)
            if gg_new != gg:
                gg       = gg_new
                remaining.discard(H)
                changed  = True
        if not changed:
            break
    return gg

def degrees(nodes, g):
    return [len(parents(g, v)) + len(children(g, v)) for v in nodes]

def sortbydegree(nodes, g):
    idx = np.argsort(degrees(nodes, g))
    return list(np.asarray(nodes)[idx])

def print_ws(ws):
    print ('{'),
    for e in ws:
        print (e, ', ')
    print ('}')

