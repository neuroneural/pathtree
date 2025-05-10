import sys, os

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np
from random import shuffle

from copy import deepcopy
from pathtree import PathTree, osumset_full
from pathtreetools import find_maximal_bcliques, full_backward_inference, to_path_forest, forest_to_set

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

def find_crd_paths(graph, source, target):
    """
    Finds all CRD paths from source to target in the graph.
    A CRD path is a directed path with no repeated cycles.

    :param graph: The input graph (dictionary format).
    :param source: The starting node.
    :param target: The ending node.
    :return: List of tuples (path, total_lag).
    """
    def dfs(current, path, lag_sum, visited_edges):
        if current == target:
            paths.append((path.copy(), lag_sum))
            return

        for neighbor in graph.get(current, {}):
            edge_lags = graph[current][neighbor].get(1, set())
            if not edge_lags:  # No valid directed edge
                continue

            # Prevent revisiting edges to break cycles
            if (current, neighbor) in visited_edges:
                continue

            for lag in edge_lags:
                path.append(neighbor)
                visited_edges.add((current, neighbor))  # Mark edge as visited
                dfs(neighbor, path, lag_sum + lag, visited_edges)
                visited_edges.remove((current, neighbor))  # Unmark edge
                path.pop()

    paths = []
    visited_edges = set()
    dfs(source, [source], 0, visited_edges)
    return paths


def find_h_treks(graph, a, b):
    """
    Finds all h-treks (heterogeneous treks) between nodes a and b in the graph.

    :param graph: The input graph (dictionary format).
    :param a: The first node (end of τ2).
    :param b: The second node (end of τ1).
    :return: List of h-treks as tuples (τ1, τ2, e).
    """
    h_treks = []

    # Find all potential roots for τ1 (paths to b)
    potential_x1 = [node for node in graph if find_crd_paths(graph, node, b)]
    # Find all potential roots for τ2 (paths to a)
    potential_x2 = [node for node in graph if find_crd_paths(graph, node, a)]

    for x1 in potential_x1:
        τ1_paths = find_crd_paths(graph, x1, b)
        for x2 in potential_x2:
            τ2_paths = find_crd_paths(graph, x2, a)

            # Case 1: τ1 and τ2 start from the same node
            if x1 == x2:
                for τ1 in τ1_paths:
                    for τ2 in τ2_paths:
                        h_treks.append((τ1, τ2, None))
            else:
                # Case 2: Check for a bi-directed edge between x1 and x2
                if x1 in graph and x2 in graph[x1] and 2 in graph[x1][x2]:
                    e = graph[x1][x2][2]
                    for τ1 in τ1_paths:
                        for τ2 in τ2_paths:
                            h_treks.append((τ1, τ2, e))
                elif x2 in graph and x1 in graph[x2] and 2 in graph[x2][x1]:
                    e = graph[x2][x1][2]
                    for τ1 in τ1_paths:
                        for τ2 in τ2_paths:
                            h_treks.append((τ1, τ2, e))

    return h_treks


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
        pt = PathTree(preset=new_base, loopset=set())
        forest = [pt]
        print("Final merged PathForest:", forest)
        return forest
    else:
        # Instead of summing just the presets, use sum_forests so that the loopset info is preserved.
        indirect_forest = sum_forests(ah_forest, hb_forest)
        # If multiple branches result, merge them into one using unify_forests.
        if len(indirect_forest) > 1:
            indirect_forest = unify_forests(indirect_forest, [])
        pt_indirect = indirect_forest[0]
        new_preset = pt_indirect.preset  # Combined preset (including any sums from ah and hb)
        
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
        
        # Combine the self-loop delays with those already in the indirect forest.
        final_loopset = pt_indirect.loopset.union(extra_hh).union(nested_hh)
        
        pt_final = PathTree(preset=new_preset, loopset=final_loopset)
        forest = [pt_final]
        print("Final merged PathForest:", forest)
        return forest
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


def merge_forest_by_base(forest):
    """
    Given a list of PathTree objects (forest), merge contributions that have the same
    (flattened) preset (assumed to be a singleton) by taking the union of their loopsets.
    
    For example, if forest contains two PathTrees both with preset {10} (and maybe differing
    loopsets), the result will be one PathTree with preset {10} and loopset equal to the union
    of the two loopsets.
    """
    merged = {}  # key: base (int), value: union of loopset contributions (a set)
    for pt in forest:
        bases = get_base_as_set(pt.preset)  # flattened preset as a plain set of ints
        loops = flatten_to_ints(pt.loopset)

        # We assume each contribution is for a single base.
        for b in bases:
            loops = flatten_to_ints(pt.loopset)
            if b in merged:
                merged[b].update(loops)
            else:
                merged[b] = set(loops)
    result = [] 
    for b, ls in merged.items():
        result.append(PathTree(preset={b}, loopset=loops))
    return result
def hide_node(g, H, zero_incoming=False):
    """
    Remove vertex H from graph g.
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
    ch = children(g, H)
    # Parents of H: dictionary {parent: delay from p->H}
    pa_orig = parents(g, H)
    if zero_incoming:
        # keep every parent, but give it zero delay
        pa = { p: {0} for p in pa_orig }
    else:
        pa = pa_orig
    # Self-loop at H:
    sl = g[H][H][1] if H in g[H] else set()

    # Remove H from gg first.
    remove_node(gg, H)

    # --- Branch 1: Update directed edges from observed parents of H to each child c ---
    if pa:
        for p in pa:
            for c in ch:
                if c in gg[p]:
                    existing = gg[p][c].get(1, None)
                else:
                    existing = None

                pa_weights = pa[p] if pa[p] else set()
                ch_weights = ch[c] if c in ch else set()
                sl_weights = sl if sl else set()

                print(f"For parent={p}, child={c}: existing={bool(existing)}, ah={pa_weights}, hb={ch_weights}, sl={sl_weights}")

                if c == H or p == H:
                    continue

                new_branch = merge_weightsets(set(), pa_weights, ch_weights, sl_weights, induced_bidirected=False)

                if existing is None:
                    gg[p][c] = {1: new_branch}
                else:
                    if not isinstance(existing, list):
                        existing = [existing]
                    combined = existing + new_branch
                    combined = merge_forest_by_base(combined)
                    gg[p][c][1] = combined

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

                # Store the induced bidirected edge as a list.
                gg[closer][further].setdefault(2, [])
                gg[closer][further][2].extend(forest)
                gg[closer][further][2] = merge_forest_by_base(gg[closer][further][2])

    # Cleanup: remove any leftover references to H in the graph.
    for node in list(gg.keys()):
        if H in gg[node]:
            del gg[node][H]
    if H in gg:
        del gg[H]

    return gg

def degrees(nodes, g):
    return [len(parents(g, v)) + len(children(g, v)) for v in nodes]


def sortbydegree(nodes, g):
    idx = np.argsort(degrees(nodes, g))
    return list(np.asarray(nodes)[idx])


def hide_nodes(g, nodelist, dosort=True):
    nodeset = set()  # make sure not to delete a node twice
    if dosort: nodelist = sortbydegree(nodelist, g)
    gg = deepcopy(g)
    for n in nodelist:
        if n in nodeset: continue
        gg = hide_node(gg, n)
        nodeset.add(n)
    return gg


def hide_random(g, ratio):
    """
    Hire random modes in the `ratio` proportion from graph g
    :param g: input graph
    :param ratio: what percentage of nodes to hide
    :return: the graph with hidden variables
    """
    nodes = g.keys()
    shuffle(nodes)
    return hide_nodes(g, nodes[:int(len(g) * ratio)])


def print_ws(ws):
    print ('{'),
    for e in ws:
        print (e, ', ')
    print ('}')

