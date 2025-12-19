"""
Test forward-backward-forward round-trip on specific graph.

Tests that:
1. Forward pass (hide nodes) with z3 implementation
2. Backward pass (reverse/expand to unit-lag with latents)
3. Second forward pass (hide the introduced latents)
4. Result should equal the first forward pass
"""

import sys
from copy import deepcopy
from pprint import pprint

from hide_nodes_z3 import hide_nodes_z3
from pathtreetools import reverse, unify_edge_representation, find_maximal_bcliques, minimal_refinement


def normalize_graph(G):
    """Normalize graph for comparison."""
    result = {}
    for u, nbrs in G.items():
        u_norm = int(u) if isinstance(u, str) and u.isdigit() else u
        if not nbrs:
            continue
        result[u_norm] = {}
        for v, ed in nbrs.items():
            v_norm = int(v) if isinstance(v, str) and v.isdigit() else v
            result[u_norm][v_norm] = {}
            for et, lags in ed.items():
                if isinstance(lags, set):
                    result[u_norm][v_norm][et] = set(lags)
                else:
                    # PathTree or other - convert to set
                    from pathtreetools import forest_to_set
                    result[u_norm][v_norm][et] = forest_to_set(lags)
    return result


def compare_lag_sets(g1, g2):
    """Compare two graphs by lag sets, return True if equal."""
    g1 = normalize_graph(g1)
    g2 = normalize_graph(g2)

    all_keys = set(g1.keys()) | set(g2.keys())

    for u in all_keys:
        nbrs1 = g1.get(u, {})
        nbrs2 = g2.get(u, {})
        all_targets = set(nbrs1.keys()) | set(nbrs2.keys())

        for v in all_targets:
            ed1 = nbrs1.get(v, {})
            ed2 = nbrs2.get(v, {})
            all_etypes = set(ed1.keys()) | set(ed2.keys())

            for et in all_etypes:
                lags1 = ed1.get(et, set())
                lags2 = ed2.get(et, set())

                if lags1 != lags2:
                    return False, (u, v, et, lags1, lags2)

    return True, None


def test_roundtrip(graph, hidden, maxlag=17, verbose=True):
    """
    Test forward-backward-forward round trip.

    Returns True if first forward == forward(backward(forward))
    """
    observed = set(graph.keys())
    for u, nbrs in graph.items():
        observed.update(nbrs.keys())
    observed -= hidden

    if verbose:
        print("=" * 60)
        print("ORIGINAL GRAPH:")
        pprint(graph)
        print(f"Hidden nodes: {hidden}")
        print(f"Observed nodes: {observed}")

    # Step 1: First forward pass (hide nodes)
    if verbose:
        print("\n" + "=" * 60)
        print("STEP 1: First forward pass (hide_nodes_z3)")

    G_forward1 = hide_nodes_z3(graph, hidden, maxlag=maxlag, verbose=False)

    if verbose:
        print("Result after first forward:")
        pprint(G_forward1)

    # Step 2: Backward pass (reverse - expand to unit-lag with latents)
    if verbose:
        print("\n" + "=" * 60)
        print("STEP 2: Backward pass (reverse)")

    # Need to unify and refine before reverse
    G_unified = unify_edge_representation(deepcopy(G_forward1))
    bcliques = find_maximal_bcliques(G_unified)
    G_refined = minimal_refinement(G_unified, bcliques)

    if verbose:
        print("Refined graph (PathTrees):")
        pprint(G_refined)

    G_backward = reverse(G_refined)

    if verbose:
        print("Result after backward (unit-lag with latents):")
        pprint(G_backward)

    # Identify introduced latents
    all_backward_nodes = set(G_backward.keys())
    for u, nbrs in G_backward.items():
        all_backward_nodes.update(nbrs.keys())

    introduced_latents = all_backward_nodes - observed

    if verbose:
        print(f"Introduced latents: {introduced_latents}")

    # Step 3: Second forward pass (hide the introduced latents)
    if verbose:
        print("\n" + "=" * 60)
        print("STEP 3: Second forward pass (hide introduced latents)")

    G_forward2 = hide_nodes_z3(G_backward, introduced_latents, maxlag=maxlag, verbose=False)

    if verbose:
        print("Result after second forward:")
        pprint(G_forward2)

    # Step 4: Compare
    if verbose:
        print("\n" + "=" * 60)
        print("COMPARISON:")

    match, diff = compare_lag_sets(G_forward1, G_forward2)

    if match:
        if verbose:
            print("✓ Round-trip successful! G_forward1 == G_forward2")
        return True
    else:
        if verbose:
            u, v, et, lags1, lags2 = diff
            etype = "directed" if et == 1 else "bidirected"
            print(f"✗ Mismatch at edge ({u}, {v}) {etype}:")
            print(f"  First forward:  {sorted(lags1)}")
            print(f"  Second forward: {sorted(lags2)}")
            print(f"  Only in first:  {sorted(lags1 - lags2)}")
            print(f"  Only in second: {sorted(lags2 - lags1)}")
        return False


def test_user_graph():
    """Test the specific graph from user request."""
    print("\n" + "=" * 70)
    print("TEST: User's graph with self-loop")
    print("=" * 70)

    # Graph with self-loop on node 5
    graph_with_loop = {
        5: {6: {1: {1}}, 7: {1: {1}}, 5: {1: {1}}},
        6: {7: {1: {1}}},
        7: {8: {1: {1}}, 6: {1: {1}}},
        8: {}
    }
    hidden = {6, 7}

    result1 = test_roundtrip(graph_with_loop, hidden, maxlag=17, verbose=True)

    print("\n" + "=" * 70)
    print("TEST: User's graph WITHOUT self-loop")
    print("=" * 70)

    # Same graph without self-loop
    graph_no_loop = {
        5: {6: {1: {1}}, 7: {1: {1}}},
        6: {7: {1: {1}}},
        7: {8: {1: {1}}, 6: {1: {1}}},
        8: {}
    }

    result2 = test_roundtrip(graph_no_loop, hidden, maxlag=17, verbose=True)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  With self-loop:    {'✓ PASS' if result1 else '✗ FAIL'}")
    print(f"  Without self-loop: {'✓ PASS' if result2 else '✗ FAIL'}")

    return result1 and result2


if __name__ == "__main__":
    success = test_user_graph()
    sys.exit(0 if success else 1)
