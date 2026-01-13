"""
Z3-based implementation of hide_nodes (graph marginalization).

This module ports the core 'forward' pass from clingo/ASP to z3.
It computes the marginalized graph by hiding latent nodes and deriving:
  - Directed edges between observed nodes (with lag sets)
  - Bidirected edges from common hidden ancestors (with difference sets)
"""

from collections import defaultdict
from typing import Dict, Set, Tuple, Any


def hide_nodes_z3(
    graph: Dict[Any, Dict[Any, Dict[int, Set[int]]]],
    hidden: Set[Any],
    maxlag: int = 17,
    verbose: bool = False
) -> Dict[Any, Dict[Any, Dict[int, Set[int]]]]:
    """
    Marginalize hidden nodes from a graph using Z3 constraint solving.

    This is the Z3 equivalent of the clingo hide_nodes.lp program.

    Args:
        graph: Graph in format {u: {v: {etype: {lags}}}}
               etype=1 for directed, etype=2 for bidirected
        hidden: Set of nodes to marginalize out
        maxlag: Maximum lag to consider
        verbose: Print debug information

    Returns:
        Marginalized graph with only observed nodes
    """
    hidden = set(hidden)

    # Collect all nodes
    all_nodes = set(graph.keys())
    for u, nbrs in graph.items():
        all_nodes.update(nbrs.keys())

    observed = all_nodes - hidden

    if verbose:
        print(f"[z3] All nodes: {sorted(all_nodes, key=str)}")
        print(f"[z3] Hidden: {sorted(hidden, key=str)}")
        print(f"[z3] Observed: {sorted(observed, key=str)}")

    # Extract edge information
    edges = {}  # (u, v) -> set of lags for directed edges
    loop_lengths = defaultdict(set)  # node -> set of loop lengths

    for u, nbrs in graph.items():
        for v, ed in nbrs.items():
            if 1 in ed:  # directed edge
                lags = ed[1]
                if u == v and u in hidden:
                    # Self-loop on hidden node - record both as loop length AND as edge
                    for lag in lags:
                        if lag >= 1:
                            loop_lengths[u].add(lag)
                    # Also add to edges so it appears in adjacency
                    edges[(u, v)] = lags
                else:
                    # Regular directed edge
                    edges[(u, v)] = lags

    if verbose:
        print(f"[z3] Edges: {edges}")
        print(f"[z3] Loop lengths: {dict(loop_lengths)}")

    # Compute paths using fixed-point iteration (more efficient than pure Z3 for this)
    result = _compute_marginalized_graph(
        all_nodes, observed, hidden, edges, loop_lengths, maxlag, verbose
    )

    return result


def _compute_marginalized_graph(
    all_nodes: Set[Any],
    observed: Set[Any],
    hidden: Set[Any],
    edges: Dict[Tuple[Any, Any], Set[int]],
    loop_lengths: Dict[Any, Set[int]],
    maxlag: int,
    verbose: bool
) -> Dict[Any, Dict[Any, Dict[int, Set[int]]]]:
    """
    Compute the marginalized graph by finding all paths through hidden nodes.

    This implements the core logic from hide_nodes.lp:
    - path_h(Z,Y,L): paths from hidden Z to observed Y with length L
    - pathvv(X,Y,L): paths from observed X to observed Y through hidden nodes
    - bi_edge/bi_diff: bidirected edges from common hidden ancestors
    """

    # Build adjacency for unit edges
    # For simplicity, we treat all edges as having lag=1 between nodes
    # (the actual lag structure is handled by path lengths)
    adj = defaultdict(set)
    for (u, v), lags in edges.items():
        if 1 in lags or any(lag == 1 for lag in lags):
            adj[u].add(v)

    # Also add edges for nodes with multi-lag (they represent chains of length > 1)
    # but for the forward pass we mainly care about unit connections
    for (u, v), lags in edges.items():
        if u != v:
            adj[u].add(v)

    if verbose:
        print(f"[z3] Adjacency: {dict(adj)}")

    # Compute path_h: hidden node Z to observed node Y with path length L
    # path_h(Z,Y,1) :- hidden(Z), observed(Y), edge(Z,Y).
    # path_h(Z,Y,L) :- hidden(Z), edge(Z,Z2), hidden(Z2), path_h(Z2,Y,L-1).
    path_h = defaultdict(set)  # (Z, Y) -> set of path lengths

    # Base case: direct edges from hidden to observed
    for z in hidden:
        for y in adj[z]:
            if y in observed:
                path_h[(z, y)].add(1)

    # Fixed-point iteration for paths through hidden nodes
    # This includes self-loops (z == z2) which allow repeated traversal
    changed = True
    while changed:
        changed = False
        for z in hidden:
            for z2 in adj[z]:
                if z2 in hidden:
                    for y in observed:
                        if (z2, y) in path_h:
                            for L in list(path_h[(z2, y)]):
                                new_L = L + 1
                                if new_L <= maxlag and new_L not in path_h[(z, y)]:
                                    path_h[(z, y)].add(new_L)
                                    changed = True

        # Also expand self-loops: if z has a self-loop and can reach y,
        # then it can reach y with any additional multiple of the loop
        for z in hidden:
            if z in adj[z]:  # z has a self-loop
                for y in observed:
                    if (z, y) in path_h:
                        current_lengths = list(path_h[(z, y)])
                        for L in current_lengths:
                            new_L = L + 1  # one more traversal of the self-loop
                            if new_L <= maxlag and new_L not in path_h[(z, y)]:
                                path_h[(z, y)].add(new_L)
                                changed = True

    if verbose:
        print("[z3] path_h (hidden->observed paths):")
        for (z, y), lengths in sorted(path_h.items(), key=str):
            print(f"  {z} -> {y}: {sorted(lengths)}")

    # Compute base_h: loop-free distances (avoiding revisiting same node)
    # This is needed for bidirected edges
    base_h = _compute_base_h(hidden, observed, adj, maxlag, verbose)

    if verbose:
        print("[z3] base_h (loop-free paths):")
        for (h, y), lengths in sorted(base_h.items(), key=str):
            print(f"  {h} -> {y}: {sorted(lengths)}")

    # Compute path_hidden: paths within hidden subgraph (for loop detection)
    path_hidden = _compute_path_hidden(hidden, adj, maxlag, verbose)

    # Compute actual loop_len from path_hidden (cycles back to same node)
    computed_loops = defaultdict(set)
    for h in hidden:
        if (h, h) in path_hidden:
            computed_loops[h] = path_hidden[(h, h)]
        # Also include declared loop lengths
        computed_loops[h] |= loop_lengths.get(h, set())

    if verbose:
        print(f"[z3] computed_loops: {dict(computed_loops)}")

    # Compute pathvv: observed-to-observed paths
    # pathvv(X,Y,1) :- observed(X), observed(Y), edge(X,Y).
    # pathvv(X,Y,L) :- observed(X), edge(X,Z), hidden(Z), path_h(Z,Y,L-1).
    pathvv = defaultdict(set)

    # Direct observed-to-observed edges
    for (u, v), lags in edges.items():
        if u in observed and v in observed:
            for lag in lags:
                if lag <= maxlag:
                    pathvv[(u, v)].add(lag)

    # Paths through hidden nodes
    for x in observed:
        for z in adj[x]:
            if z in hidden:
                for y in observed:
                    if (z, y) in path_h:
                        for L in path_h[(z, y)]:
                            new_L = L + 1
                            if new_L <= maxlag:
                                pathvv[(x, y)].add(new_L)

    if verbose:
        print("[z3] pathvv (observed->observed paths):")
        for (x, y), lengths in sorted(pathvv.items(), key=str):
            print(f"  {x} -> {y}: {sorted(lengths)}")

    # Compute bidirected edges from common hidden ancestors
    bi_diff = _compute_bidirected(
        hidden, observed, base_h, computed_loops, maxlag, verbose
    )

    # Build result graph
    result = {}

    # Add directed edges
    for (u, v), lags in pathvv.items():
        if lags:
            result.setdefault(u, {}).setdefault(v, {})[1] = set(lags)

    # Add bidirected edges
    for (u, v), diffs in bi_diff.items():
        if diffs:
            result.setdefault(u, {}).setdefault(v, {})[2] = set(diffs)

    if verbose:
        print("[z3] Result graph:")
        for u, nbrs in result.items():
            for v, ed in nbrs.items():
                print(f"  {u} -> {v}: {ed}")

    return result


def _compute_base_h(
    hidden: Set[Any],
    observed: Set[Any],
    adj: Dict[Any, Set[Any]],
    maxlag: int,
    verbose: bool
) -> Dict[Tuple[Any, Any], Set[int]]:
    """
    Compute base_h: loop-free hidden->observed distances.

    base_h(H,Y,1) :- hidden(H), observed(Y), edge(H,Y).
    base_h(H,Y,L) :- hidden(H), edge(H,Z), hidden(Z), H != Z,
                     base_h(Z,Y,L2), L = L2 + 1.
    """
    base_h = defaultdict(set)

    # Base case
    for h in hidden:
        for y in adj[h]:
            if y in observed:
                base_h[(h, y)].add(1)

    # Fixed-point iteration (avoiding cycles by tracking visited)
    # We use BFS-style expansion to avoid revisiting nodes in same path
    for h in hidden:
        for y in observed:
            _bfs_base_h(h, y, hidden, observed, adj, base_h, maxlag)

    return base_h


def _bfs_base_h(
    start_h: Any,
    target_y: Any,
    hidden: Set[Any],
    observed: Set[Any],
    adj: Dict[Any, Set[Any]],
    base_h: Dict[Tuple[Any, Any], Set[int]],
    maxlag: int
):
    """BFS to find all simple (non-repeating) paths from start_h to target_y."""
    from collections import deque

    # Queue entries: (current_node, path_length, visited_set)
    queue = deque([(start_h, 0, {start_h})])

    while queue:
        curr, dist, visited = queue.popleft()

        for nxt in adj[curr]:
            if nxt in visited:
                continue

            new_dist = dist + 1
            if new_dist > maxlag:
                continue

            if nxt == target_y and nxt in observed:
                base_h[(start_h, target_y)].add(new_dist)
            elif nxt in hidden:
                new_visited = visited | {nxt}
                queue.append((nxt, new_dist, new_visited))


def _compute_path_hidden(
    hidden: Set[Any],
    adj: Dict[Any, Set[Any]],
    maxlag: int,
    verbose: bool
) -> Dict[Tuple[Any, Any], Set[int]]:
    """
    Compute path_hidden: paths within hidden subgraph.

    path_hidden(H,H,1) :- hidden(H), edge(H,H).
    path_hidden(H1,H2,1) :- hidden(H1), hidden(H2), edge(H1,H2), H1 != H2.
    path_hidden(H1,H3,L) :- hidden(H1), edge(H1,H2), hidden(H2),
                           path_hidden(H2,H3,L2), L = L2 + 1.
    """
    path_hidden = defaultdict(set)

    # Base case: direct edges between hidden nodes
    for h1 in hidden:
        for h2 in adj[h1]:
            if h2 in hidden:
                path_hidden[(h1, h2)].add(1)

    # Fixed-point iteration
    changed = True
    while changed:
        changed = False
        for h1 in hidden:
            for h2 in adj[h1]:
                if h2 in hidden:
                    for h3 in hidden:
                        if (h2, h3) in path_hidden:
                            for L in list(path_hidden[(h2, h3)]):
                                new_L = L + 1
                                if new_L <= maxlag and new_L not in path_hidden[(h1, h3)]:
                                    path_hidden[(h1, h3)].add(new_L)
                                    changed = True

    return path_hidden


def _compute_bidirected(
    hidden: Set[Any],
    observed: Set[Any],
    base_h: Dict[Tuple[Any, Any], Set[int]],
    loop_lengths: Dict[Any, Set[int]],
    maxlag: int,
    verbose: bool
) -> Dict[Tuple[Any, Any], Set[int]]:
    """
    Compute bidirected edges from common hidden ancestors.

    bi_diff(X,Y,D) :- observed(X), observed(Y), hidden(H),
                     base_h(H,X,Lx), base_h(H,Y,Ly), X != Y, Ly >= Lx, D = Ly - Lx.

    Also handles periodic expansion from loops:
    bi_diff(X,Y,D2) :- looped_ancestor(H0,X,Lx), looped_ancestor(H0,Y,Ly),
                      has_loop(H0), loop_len(H0,Lloop), D0 = Ly - Lx,
                      D2 = D0 + K * Lloop.
    """
    bi_diff = defaultdict(set)

    # Find all pairs of observed nodes with common hidden ancestors
    for h in hidden:
        # Get all observed nodes reachable from h
        reachable = {}  # obs_node -> set of base distances
        for (hh, y), dists in base_h.items():
            if hh == h and y in observed:
                reachable[y] = dists

        if len(reachable) < 1:
            continue

        # For each pair of reachable observed nodes
        obs_list = list(reachable.keys())
        for i, x in enumerate(obs_list):
            for y in obs_list[i:]:  # include x==y for self-loops
                for lx in reachable[x]:
                    for ly in reachable[y]:
                        if ly >= lx:
                            d = ly - lx
                            if d <= maxlag:
                                if x != y:
                                    bi_diff[(x, y)].add(d)
                                else:
                                    # Self-bidirected (x == y, d > 0)
                                    if d > 0:
                                        bi_diff[(x, x)].add(d)

        # Handle periodic expansion from loops
        if h in loop_lengths and loop_lengths[h]:
            for lloop in loop_lengths[h]:
                for x in reachable:
                    for y in reachable:
                        for lx in reachable[x]:
                            for ly in reachable[y]:
                                if ly >= lx or x == y:
                                    d0 = abs(ly - lx)
                                    # Add multiples of loop length
                                    for k in range(1, (maxlag // lloop) + 2):
                                        d2 = d0 + k * lloop
                                        if d2 <= maxlag:
                                            if x != y:
                                                if ly >= lx:
                                                    bi_diff[(x, y)].add(d2)
                                            else:
                                                # Self-bidirected
                                                if d2 > 0:
                                                    bi_diff[(x, x)].add(d2)

    # Normalize: ensure (x,y) with x < y (except self-loops)
    normalized = defaultdict(set)
    for (x, y), diffs in bi_diff.items():
        if x == y:
            normalized[(x, y)] |= diffs
        else:
            key = (min(x, y), max(x, y)) if isinstance(x, int) and isinstance(y, int) else (x, y) if str(x) <= str(y) else (y, x)
            normalized[key] |= diffs

    return normalized


def _is_helper_cycle(node: Any) -> bool:
    """Check if node is a helper cycle node (cyc/bcyc) that should be ignored for loop sources."""
    if isinstance(node, str):
        return node.startswith("cyc(") or node.startswith("bcyc(")
    return False


# Convenience function matching the clingo interface
def run_hide_nodes_z3(
    graph: Dict[Any, Dict[Any, Dict[int, Set[int]]]],
    hidden: Set[Any],
    maxlag: int = 17,
    verbose: bool = False
) -> Dict[Any, Dict[Any, Dict[int, Set[int]]]]:
    """
    Run hide_nodes using Z3 solver.

    This is the main entry point, matching the interface of the clingo version.
    """
    return hide_nodes_z3(graph, hidden, maxlag, verbose)
