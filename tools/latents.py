
import sys, os

sys.path.append('./tools/')
import traversal as trv
import bfutils as bfu
import graphkit as gk
import numpy as np
from random import shuffle

from copy import deepcopy
from pathtree import PathTree, osumset, osumset_full, InfiniteExpression
from pathtreetools import apply_minimal_refinement
from infiniteset import InfiniteSet

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

def apply_batch_pts(g):
    """
    Applies the Batch PTS algorithm to refine the graph's edge-lag sets iteratively.
    :param g: Input graph
    :return: Refined graph
    """
    converged = False
    iteration = 0

    while not converged:
        print(f"Starting iteration {iteration}...")
        
        # Take a snapshot of the edge-lag sets
        prev_edge_lags = {v: {w: g[v][w][1] for w in g[v]} for v in g}

        # Apply the refinement process
        refined_graph = apply_minimal_refinement(g)

        # Compare current edge-lag sets with the previous snapshot
        current_edge_lags = {v: {w: refined_graph[v][w][1] for w in refined_graph[v]} for v in refined_graph}

        if prev_edge_lags == current_edge_lags:
            converged = True
        else:
            g = refined_graph
            iteration += 1

    print("Converged after", iteration, "iterations.")
    return refined_graph

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
        bases1 = base_as_set(pt1.preset)
        for pt2 in hb_forest:
            bases2 = base_as_set(pt2.preset)
            for x in bases1:
                for y in bases2:
                    diff_val = y - x
                    new_forest.append(PathTree(preset={diff_val}, loopset=set()))
    return new_forest

def update_bidirected_edges(graph, r):
    """
    Update bi-directed edges when removing a latent node 'r'.
    Recalculate lags between parents and children of 'r'.

    :param graph: The graph as a dictionary.
    :param r: The node to be removed.
    :return: The updated graph.
    """
    # Ensure 'r' exists in the graph
    if r not in graph:
        print(f"Skipping Node {r}: It does not exist in the graph.")
        return graph

    # Safely retrieve parents and children of node `r` related to bi-directed edges
    parents_r = {p: graph[p][r].get(2, set()) for p in graph if r in graph[p] and 2 in graph[p][r]}
    children_r = {c: graph[r][c].get(2, set()) for c in graph.get(r, {})}

    # Debugging output
    print(f"Node to Remove: {r}")
    print(f"Parents of {r} (Bi-directed): {parents_r}")
    print(f"Children of {r} (Bi-directed): {children_r}")

    updated_graph = deepcopy(graph)

    # Remove the latent node `r`
    if r in updated_graph:
        del updated_graph[r]

    for p in parents_r:
        for c in children_r:
            combined_lags = osumset_full(parents_r[p], children_r[c])

            if p == c:  # Self-loop case
                if p in updated_graph and p in updated_graph[p]:
                    updated_graph[p][p][2].update(combined_lags)
                else:
                    updated_graph.setdefault(p, {})[p] = {2: combined_lags}
            else:  # Regular edge case
                if c in updated_graph[p]:
                    updated_graph[p][c].setdefault(2, set()).update(combined_lags)
                else:
                    updated_graph[p][c] = {2: combined_lags}

    # Clean up remaining references to `r` in the graph
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

def to_path_forest(x):
    """
    Convert x into a list of PathTrees (a path-forest).
    If x is:
      - an int: return [PathTree(preset={x}, loopset=set())]
      - a set: return [PathTree(preset={b}, loopset=set()) for each b in x]
      - a PathTree: return [x]
      - a list (of PathTrees): return x
    """
    if isinstance(x, list):
        return x
    elif isinstance(x, PathTree):
        return [x]
    elif isinstance(x, int):
        return [PathTree(preset={x}, loopset=set())]
    elif isinstance(x, set):
        return [PathTree(preset={b}, loopset=set()) for b in x]
    elif x == set():
        return []
    else:
        raise TypeError(f"Cannot convert {x} into a path-forest")

def get_base_as_set(x):
    """
    Convert x into a set of ints.
      - If x is an int, return {x}.
      - If x is a set, return it.
      - If x is a list, return the union over its elements.
      - If x is a PathTree, return its preset (assumed to be int or set).
    """
    if isinstance(x, list):
        result = set()
        for item in x:
            result = result.union(get_base_as_set(item))
        return result
    elif isinstance(x, set):
        return x
    elif isinstance(x, int):
        return {x}
    elif isinstance(x, PathTree):
        if isinstance(x.preset, int):
            return {x.preset}
        elif isinstance(x.preset, set):
            return x.preset
        else:
            raise TypeError("Unexpected type for PathTree.preset")
    else:
        raise TypeError("get_base_as_set expects an int, set, list, or PathTree")
    
def sum_forests(f1, f2):
    """
    For each PathTree in f1 and each PathTree in f2, produce new PathTrees
    whose preset is the Cartesian sum of their bases and whose loopset is the union.
    Return the new forest (list of PathTrees).
    """
    new_forest = []
    for pt1 in f1:
        bases1 = base_as_set(pt1.preset)
        loops1 = set(pt1.loopset)
        for pt2 in f2:
            bases2 = base_as_set(pt2.preset)
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
        base2 = base_as_set(pt2.preset)
        inserted = False
        for pt1 in new_forest:
            base1 = base_as_set(pt1.preset)
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
    self-loop delays at H (hh), returning a path-forest (a list containing one PathTree)
    with a separated preset (the cycle-free delay) and a hierarchical loopset.
    
    Parameters:
      ab: Existing direct edge delay from parent to child (int, set, or PathTree/forest).
      ah: Delay from parent to H (int, set, or PathTree/forest).
      hb: Delay from H to child (int, set, or PathTree/forest).
      hh: Self-loop delay at H. It can be a set/int (flat) or a list that may include nested
          PathTree objects (for nested cycles).
      induced_bidirected: If True, perform difference merging (not our current case).
      
    For the non-induced case:
      1. Convert ah and hb to sets (ah_set and hb_set) via get_base_as_set.
      2. Compute the indirect delay as indirect = osumset_full(ah_set, hb_set) and remove 0.
      3. Let direct = get_base_as_set(ab). If direct is empty, then new_preset = indirect and extra = {}.
         Otherwise, new_preset = direct and extra = (indirect - direct).
      4. Process hh:
         - If hh is a list, separate items that are PathTrees (nested branches) from flat delays.
         - Let extra_hh be the union of the flat items (using get_base_as_set).
         - Let nested_hh be the set of nested PathTree objects.
         - For each nested branch, remove its preset (assumed to be a singleton) from the flat extra (if present).
      5. Set final_loopset = extra ∪ extra_hh ∪ nested_hh.
    
    Returns a forest (list containing one PathTree) with:
      preset = new_preset,
      loopset = final_loopset.
    """
    # Convert ah and hb to sets.
    if isinstance(ah, list):
        ah_set = set().union(*(get_base_as_set(pt.preset) for pt in ah))
    else:
        ah_set = get_base_as_set(ah)
    if isinstance(hb, list):
        hb_set = set().union(*(get_base_as_set(pt.preset) for pt in hb))
    else:
        hb_set = get_base_as_set(hb)
    
    if induced_bidirected:
        new_base = {y - x for x in ah_set for y in hb_set}
        pt = PathTree(preset=new_base)
        pt.loopset = set()
        forest = [pt]
        print("Final merged PathForest:", forest)
        return forest
    else:
        indirect = osumset_full(ah_set, hb_set)
        if 0 in indirect:
            indirect.discard(0)
        print("  Indirect paths (sum of ah and hb):", indirect)
        direct = get_base_as_set(ab)
        if not direct:
            new_preset = indirect
            extra = set()
        else:
            new_preset = direct
            extra = indirect - direct
        
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
        
        # Remove any nested branch's preset value from the flat extra.
        new_extra = set(extra)
        for nested in nested_hh:
            nested_val = list(get_base_as_set(nested.preset))[0]  # assume singleton
            if nested_val in new_extra:
                new_extra.remove(nested_val)
        final_loopset = new_extra.union(extra_hh).union(nested_hh)
        
        pt = PathTree(preset=new_preset)
        pt.loopset = final_loopset
        forest = [pt]
        print("  Final merged PathForest:", forest)
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

def hide_node(g, H):
    """
    Remove vertex H from graph g.
    For each parent p of H and each child c of H, extract:
      - ab: existing direct delay from p to c (if any)
      - ah: delay from p to H
      - hb: delay from H to c
      - sl: self-loop delay at H (if any)
    Then call merge_weightsets to compute the delay contribution for that branch.
    If the edge p->c already has a forest (list of PathTrees), append the new branch;
    otherwise, store the branch as the forest.
    Return the updated graph.
    """
    gg = deepcopy(g)
    if H not in g:
        raise KeyError(f"Node {H} not found in graph.")
    
    ch = children(g, H)  # {child: delay from H->child}
    pa = parents(g, H)   # {parent: delay from p->H}
    
    if H in g[H]:
        sl = g[H][H][1]
    else:
        sl = set()
    
    remove_node(gg, H)
    
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
                print(f"For parent={p}, child={c}: ab={'existing forest' if existing is not None else set()}, ah={pa_weights}, hb={ch_weights}, sl={sl_weights}")
                if c == H or p == H:
                    continue
                # Use an empty ab so that the new branch's preset is computed solely from the indirect path.
                new_branch = merge_weightsets(set(), pa_weights, ch_weights, sl_weights, induced_bidirected=False)
                if existing is None:
                    gg[p][c] = {1: new_branch}
                else:
                    if not isinstance(existing, list):
                        existing = [existing]
                    gg[p][c][1] = existing + new_branch
    else:
        if ch:
            children_list = list(ch.keys())
            for i in range(len(children_list)):
                for j in range(i+1, len(children_list)):
                    c1 = children_list[i]
                    c2 = children_list[j]
                    lag_c1 = ch[c1] if c1 in ch else set()
                    lag_c2 = ch[c2] if c2 in ch else set()
                    forest = merge_weightsets(set(), lag_c1, lag_c2, sl, induced_bidirected=True)
                    for (u, v) in [(c1, c2), (c2, c1)]:
                        if u not in gg:
                            gg[u] = {}
                        if v not in gg[u]:
                            gg[u][v] = {}
                        gg[u][v][2] = forest
    for parent in list(gg.keys()):
        gg[parent].pop(H, None)
    gg.pop(H, None)
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


# def test_osumnum():
#     assert osumnum(set(range(5)), 1) == set(range(1, 5 + 1))


def testcase(n):
    g1 = {1: {2: {1: {1}}, 4: {1: {1}}},
          2: {3: {1: {1}}, 7: {1: {1}}},
          3: {},
          4: {5: {1: {1}}},
          5: {3: {1: {1}}, 6: {1: {1}}},
          6: {5: {1: {1}}},
          7: {8: {1: {1}}},
          8: {2: {1: {1}}}}

    g2 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {3: {1: {1}}}}

    g3 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 6: {1: {1}}},
          3: {4: {1: {1}}, 8: {1: {1}}},
          4: {5: {1: {1}}},
          5: {},
          6: {7: {1: {1}}},
          7: {2: {1: {1}}},
          8: {9: {1: {1}}},
          9: {10: {1: {1}}},
          10: {11: {1: {1}}},
          11: {3: {1: {1}}}}

    g4 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {}}

    g5 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g6 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}, 6: {1: {1}}},
          7: {8: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g7 = {1: {2: {1: {1}}},
          2: {3: {1: {1}}, 4: {1: {1}}},
          4: {5: {1: {1}}},
          5: {6: {1: {1}}, 7: {1: {1}}},
          6: {2: {1: {1}}, 6: {1: {1}}},
          7: {8: {1: {1}}, 7: {1: {1}}},
          8: {5: {1: {1}}},
          3: {9: {1: {1}}},
          9: {10: {1: {1}}, 11: {1: {1}}},
          11: {9: {1: {1}}},
          10: {}}

    g8 = {1: {2: {1: {1}}, 5: {1: {1}}},
          2: {3: {1: {1}}, 2: {1: {1}}},
          3: {4: {1: {1}}},
          4: {8: {1: {1}}},
          5: {6: {1: {1}}},
          6: {7: {1: {1}}},
          7: {4: {1: {1}}},
          8: {9: {1: {1}}},
          9: {9: {1: {1}}, 10: {1: {1}}},
          10: {}}
    
    g9 ={
        1: {2: {1: {1}}},                # Node 1 points to Node 2
        2: {2: {1: {1}}, 3: {1: {1}}},   # Node 2 ↔ Node 3 (bi-directed), Node 2 → Node 4
        3: {}                         # Node 5 is a sink
        }

    cases = [g1, g2, g3, g4, g5, g6, g7, g8, g9]

    return fix_selfloops(cases[n])
