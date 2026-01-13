"""
Test script to verify Z3 hide_nodes produces correct outputs.

Since clingo may not be available, we test against known expected results.
"""

import sys
from pprint import pprint

# Import Z3 implementation
from hide_nodes_z3 import hide_nodes_z3


def normalize_graph(G):
    """Normalize graph for comparison (sort sets, normalize keys)."""
    result = {}
    for u, nbrs in G.items():
        u_norm = int(u) if isinstance(u, str) and u.isdigit() else u
        result[u_norm] = {}
        for v, ed in nbrs.items():
            v_norm = int(v) if isinstance(v, str) and v.isdigit() else v
            result[u_norm][v_norm] = {}
            for et, lags in ed.items():
                result[u_norm][v_norm][et] = set(lags)
    return result


def compare_graphs(g1, g2, name1="Expected", name2="Got"):
    """Compare two graphs and report differences."""
    g1 = normalize_graph(g1)
    g2 = normalize_graph(g2)

    all_keys = set(g1.keys()) | set(g2.keys())
    differences = []

    for u in sorted(all_keys, key=str):
        nbrs1 = g1.get(u, {})
        nbrs2 = g2.get(u, {})
        all_targets = set(nbrs1.keys()) | set(nbrs2.keys())

        for v in sorted(all_targets, key=str):
            ed1 = nbrs1.get(v, {})
            ed2 = nbrs2.get(v, {})
            all_etypes = set(ed1.keys()) | set(ed2.keys())

            for et in sorted(all_etypes):
                lags1 = ed1.get(et, set())
                lags2 = ed2.get(et, set())

                if lags1 != lags2:
                    etype_name = "directed" if et == 1 else "bidirected"
                    differences.append({
                        "edge": (u, v),
                        "etype": etype_name,
                        f"{name1}": sorted(lags1),
                        f"{name2}": sorted(lags2),
                        "only_in_expected": sorted(lags1 - lags2),
                        "only_in_got": sorted(lags2 - lags1),
                    })

    return differences


def test_simple_chain():
    """Test: simple chain 1 -> 2 -> 3 with node 2 hidden."""
    print("\n" + "=" * 60)
    print("TEST: Simple chain (1 -> 2 -> 3, hide {2})")
    print("=" * 60)

    # Graph: 1 -> 2 -> 3 (all unit lag)
    graph = {
        1: {2: {1: {1}}},
        2: {3: {1: {1}}},
    }
    hidden = {2}

    # Expected: 1 -> 3 with lag 2 (path through hidden node 2)
    expected = {
        1: {3: {1: {2}}}
    }

    g_z3 = hide_nodes_z3(graph, hidden=hidden, maxlag=10, verbose=False)

    print("\nExpected:")
    pprint(expected)
    print("\nZ3 result:")
    pprint(g_z3)

    diffs = compare_graphs(expected, g_z3)
    if diffs:
        print("\nDifferences found:")
        for d in diffs:
            pprint(d)
        return False
    else:
        print("\n✓ Results match!")
        return True


def test_diamond():
    """Test: diamond pattern with hidden middle nodes."""
    print("\n" + "=" * 60)
    print("TEST: Diamond (1 -> {2,3} -> 4, hide {2,3})")
    print("=" * 60)

    # Graph: 1 -> 2 -> 4, 1 -> 3 -> 4
    graph = {
        1: {2: {1: {1}}, 3: {1: {1}}},
        2: {4: {1: {1}}},
        3: {4: {1: {1}}},
    }
    hidden = {2, 3}

    # Expected: 1 -> 4 with lag 2 (two paths, both length 2)
    # Also bidirected 2 <-> 3 from common ancestor 1? No, 1 is observed.
    # But 2 and 3 are hidden, so no bidirected for them.
    expected = {
        1: {4: {1: {2}}}
    }

    g_z3 = hide_nodes_z3(graph, hidden=hidden, maxlag=10, verbose=False)

    print("\nExpected:")
    pprint(expected)
    print("\nZ3 result:")
    pprint(g_z3)

    diffs = compare_graphs(expected, g_z3)
    if diffs:
        print("\nDifferences found:")
        for d in diffs:
            pprint(d)
        return False
    else:
        print("\n✓ Results match!")
        return True


def test_common_ancestor():
    """Test: common hidden ancestor creating bidirected edge."""
    print("\n" + "=" * 60)
    print("TEST: Common ancestor (3 -> 1, 3 -> 2, hide {3})")
    print("=" * 60)

    # Graph: hidden node 3 points to both 1 and 2
    graph = {
        3: {1: {1: {1}}, 2: {1: {1}}},
    }
    hidden = {3}

    # Expected: bidirected edge 1 <-> 2 with diff=0 (same distance from 3)
    expected = {
        1: {2: {2: {0}}}
    }

    g_z3 = hide_nodes_z3(graph, hidden=hidden, maxlag=10, verbose=False)

    print("\nExpected:")
    pprint(expected)
    print("\nZ3 result:")
    pprint(g_z3)

    diffs = compare_graphs(expected, g_z3)
    if diffs:
        print("\nDifferences found:")
        for d in diffs:
            pprint(d)
        return False
    else:
        print("\n✓ Results match!")
        return True


def test_loop():
    """Test: hidden node with self-loop."""
    print("\n" + "=" * 60)
    print("TEST: Loop (1 -> 2 -> 2 -> 3, hide {2})")
    print("=" * 60)

    # Graph: 1 -> 2 -> 3, with 2 having a self-loop
    graph = {
        1: {2: {1: {1}}},
        2: {2: {1: {1}}, 3: {1: {1}}},
    }
    hidden = {2}

    # Expected: 1 -> 3 with lags {2, 3, 4, 5, ...} up to maxlag
    # Because path can be: 1->2->3 (lag 2), 1->2->2->3 (lag 3), etc.
    # Also: bidirected self-loop 3<->3 because hidden node 2 reaches 3
    # via multiple different path lengths, creating differences.
    expected_lags = set(range(2, 11))  # maxlag=10
    bidir_lags = set(range(1, 11))  # differences: |1-2|=1, |1-3|=2, etc.
    expected = {
        1: {3: {1: expected_lags}},
        3: {3: {2: bidir_lags}}
    }

    g_z3 = hide_nodes_z3(graph, hidden=hidden, maxlag=10, verbose=False)

    print("\nExpected:")
    pprint(expected)
    print("\nZ3 result:")
    pprint(g_z3)

    diffs = compare_graphs(expected, g_z3)
    if diffs:
        print("\nDifferences found:")
        for d in diffs:
            pprint(d)
        return False
    else:
        print("\n✓ Results match!")
        return True


def test_ringmore_small():
    """Test with a small ring graph."""
    print("\n" + "=" * 60)
    print("TEST: Small ring (5 nodes, hide {3})")
    print("=" * 60)

    # Ring: 1 -> 2 -> 3 -> 4 -> 5 -> 1
    graph = {
        1: {2: {1: {1}}},
        2: {3: {1: {1}}},
        3: {4: {1: {1}}},
        4: {5: {1: {1}}},
        5: {1: {1: {1}}},
    }
    hidden = {3}

    # Expected:
    # - Direct edges preserved: 1->2, 4->5, 5->1
    # - New edge: 2->4 with lag 2 (path through hidden 3)
    expected = {
        1: {2: {1: {1}}},
        2: {4: {1: {2}}},
        4: {5: {1: {1}}},
        5: {1: {1: {1}}},
    }

    g_z3 = hide_nodes_z3(graph, hidden=hidden, maxlag=10, verbose=False)

    print("\nExpected:")
    pprint(expected)
    print("\nZ3 result:")
    pprint(g_z3)

    diffs = compare_graphs(expected, g_z3)
    if diffs:
        print("\nDifferences found:")
        for d in diffs:
            pprint(d)
        return False
    else:
        print("\n✓ Results match!")
        return True


def test_asymmetric_ancestor():
    """Test: hidden ancestor with asymmetric distances."""
    print("\n" + "=" * 60)
    print("TEST: Asymmetric ancestor (3 -> 4 -> 1, 3 -> 2, hide {3,4})")
    print("=" * 60)

    # Graph: 3 -> 4 -> 1, 3 -> 2
    # So 3 is at distance 2 from 1, distance 1 from 2
    graph = {
        3: {4: {1: {1}}, 2: {1: {1}}},
        4: {1: {1: {1}}},
    }
    hidden = {3, 4}

    # Expected: bidirected 1 <-> 2 with diff=1 (distances 2 and 1 from ancestor 3)
    expected = {
        1: {2: {2: {1}}}
    }

    g_z3 = hide_nodes_z3(graph, hidden=hidden, maxlag=10, verbose=False)

    print("\nExpected:")
    pprint(expected)
    print("\nZ3 result:")
    pprint(g_z3)

    diffs = compare_graphs(expected, g_z3)
    if diffs:
        print("\nDifferences found:")
        for d in diffs:
            pprint(d)
        return False
    else:
        print("\n✓ Results match!")
        return True


def test_self_bidirected():
    """Test: self-bidirected edge from hidden ancestor."""
    print("\n" + "=" * 60)
    print("TEST: Self-bidirected (H -> X via two paths)")
    print("=" * 60)

    # Graph: Hidden H has two paths to observed X with different lengths
    # H -> A -> X (length 2) and H -> X (length 1)
    # This creates a self-bidirected edge X <-> X with diff=1
    graph = {
        'H': {'A': {1: {1}}, 'X': {1: {1}}},
        'A': {'X': {1: {1}}},
    }
    hidden = {'H', 'A'}

    # Expected: self-bidirected X <-> X with diff=1 (distances 1 and 2 from H)
    expected = {
        'X': {'X': {2: {1}}}
    }

    g_z3 = hide_nodes_z3(graph, hidden=hidden, maxlag=10, verbose=False)

    print("\nExpected:")
    pprint(expected)
    print("\nZ3 result:")
    pprint(g_z3)

    diffs = compare_graphs(expected, g_z3)
    if diffs:
        print("\nDifferences found:")
        for d in diffs:
            pprint(d)
        return False
    else:
        print("\n✓ Results match!")
        return True


def run_all_tests():
    """Run all tests and report summary."""
    tests = [
        test_simple_chain,
        test_diamond,
        test_common_ancestor,
        test_loop,
        test_ringmore_small,
        test_asymmetric_ancestor,
        test_self_bidirected,
    ]

    results = []
    for test in tests:
        try:
            result = test()
            results.append((test.__name__, result))
        except Exception as e:
            print(f"\n✗ {test.__name__} FAILED with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test.__name__, False))

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    passed = sum(1 for _, r in results if r)
    total = len(results)

    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {status}: {name}")

    print(f"\nTotal: {passed}/{total} tests passed")

    return passed == total


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
