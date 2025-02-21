
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
    from copy import deepcopy

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

def merge_weightsets(ab, ah, hb, hh):
    print(f"merge_weightsets called with:\n  ab={ab}, ah={ah}, hb={hb}, hh={hh}")
    
    # If there is no pre-existing direct edge (ab is empty), we are inducing a bi-directed edge.
    if not ab:
        # Compute the difference: for every parent's lag (x) and child's lag (y), compute (y - x).
        diff = {y - x for x in ah for y in hb}
        final_result = diff
    else:
        # Otherwise, use the existing summation approach.
        indirect_paths = osumset_full(ah, hb)
        print("  Indirect paths (sum of ah and hb):", indirect_paths)
        indirect_paths.discard(0)
        print("  Indirect paths after discarding 0:", indirect_paths)
        self_loop_expansion = InfiniteExpression(0)
        for h in hh:
            var = get_fresh_loop_var()
            self_loop_expansion.add_term(h, var)
        print("  Self-loop expansion:", self_loop_expansion)
        indirect_paths_with_loops = osumset_full(indirect_paths, self_loop_expansion)
        final_result = indirect_paths_with_loops.union(ab)
    
    #final_result.discard(0)
    print("  Final merged edge-lag set (ws_final):", final_result)
    return final_result if final_result else {0}


def hide_node(g, H):
    """
    Removes a node H from graph g, recalculating edges (and weights) 
    while ensuring logical handling of loops and paths.
    If H has no parents but does have children, induce new edges among its children.
    
    :param g: input graph (dictionary)
    :param H: node to hide
    :return: the updated graph with H removed
    """
    from copy import deepcopy
    gg = deepcopy(g)
    
    if H not in g:
        raise KeyError(f"Node {H} not found in graph.")
    
    # Get children and parents of the node to be hidden
    ch = children(g, H)  # e.g. { child_node: child_lags }
    pa = parents(g, H)   # e.g. { parent_node: parent_lags }
    
    # Check if H has a self-loop
    if H in g[H]:
        sl = g[H][H][1]  # e.g. {10} for a self-loop of length 10
    else:
        sl = set()
    
    # Remove the node H from gg
    remove_node(gg, H)   # deletes H and all its incident edges
    
    # If H has parents, do the usual merging.
    if pa:
        for p in pa:
            for c in ch:
                # 1) Existing p->c edge-lag set (if any)
                ab = gg[p][c].get(1, set()) if c in gg[p] else set()
                # 2) parent's p->H lag set
                pa_weights = pa[p] if pa[p] else set()
                # 3) H->child c lag set
                ch_weights = ch[c] if c in ch else set()
                # 4) self-loop on H
                sl_weights = sl if sl else set()
                print(f"For parent={p}, child={c}: ab={ab}, ah={pa_weights}, hb={ch_weights}, sl={sl_weights}")
                if c == H or p == H:
                    continue
                w = merge_weightsets(ab, pa_weights, ch_weights, sl_weights)
                if p not in gg:
                    gg[p] = {}
                if c not in gg[p]:
                    gg[p][c] = {}
                gg[p][c][1] = w
    else:
        # If H has no parents (i.e., H is a source) but has children,
        # then for every pair of distinct children, create an induced edge.
        if ch:
            children_list = list(ch.keys())
            for i in range(len(children_list)):
                for j in range(i+1, len(children_list)):
                    c1 = children_list[i]
                    c2 = children_list[j]
                    lag_c1 = ch[c1] if c1 in ch else set()
                    lag_c2 = ch[c2] if c2 in ch else set()
                    w = merge_weightsets(set(), lag_c1, lag_c2, sl)
                    # Insert the induced edge as a bi-directed edge (edge type 2) in both directions.
                    for (u, v) in [(c1, c2), (c2, c1)]:
                        if u not in gg:
                            gg[u] = {}
                        if v not in gg[u]:
                            gg[u][v] = {}
                        gg[u][v][2] = w
    # Clean up any remaining references to H.
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
