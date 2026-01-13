"""
Pytest tests for hide_nodes_z3 (forward pass).

Tests graph marginalization with various topologies:
- Simple chains
- Diamond patterns
- Common ancestors (bidirected edges)
- Self-loops on hidden nodes
- Ring structures
- Asymmetric distances
- Self-bidirected edges
"""

from hide_nodes_z3 import hide_nodes_z3


def normalize_graph(G):
    """Normalize graph for comparison."""
    result = {}
    for u, nbrs in G.items():
        u_norm = int(u) if isinstance(u, str) and str(u).isdigit() else u
        result[u_norm] = {}
        for v, ed in nbrs.items():
            v_norm = int(v) if isinstance(v, str) and str(v).isdigit() else v
            result[u_norm][v_norm] = {}
            for et, lags in ed.items():
                result[u_norm][v_norm][et] = set(lags)
    return result


def graphs_equal(g1, g2):
    """Check if two graphs have equal edge sets."""
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
                if ed1.get(et, set()) != ed2.get(et, set()):
                    return False
    return True


class TestSimpleChains:
    """Test simple chain topologies."""

    def test_two_node_chain(self):
        """1 -> 2 -> 3, hide {2} => 1 -> 3 with lag 2"""
        graph = {
            1: {2: {1: {1}}},
            2: {3: {1: {1}}},
        }
        hidden = {2}
        expected = {1: {3: {1: {2}}}}

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, expected)

    def test_three_node_chain(self):
        """1 -> 2 -> 3 -> 4, hide {2, 3} => 1 -> 4 with lag 3"""
        graph = {
            1: {2: {1: {1}}},
            2: {3: {1: {1}}},
            3: {4: {1: {1}}},
        }
        hidden = {2, 3}
        expected = {1: {4: {1: {3}}}}

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, expected)


class TestDiamondPatterns:
    """Test diamond/fork topologies."""

    def test_diamond(self):
        """1 -> {2,3} -> 4, hide {2,3} => 1 -> 4 with lag 2"""
        graph = {
            1: {2: {1: {1}}, 3: {1: {1}}},
            2: {4: {1: {1}}},
            3: {4: {1: {1}}},
        }
        hidden = {2, 3}
        expected = {1: {4: {1: {2}}}}

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, expected)


class TestBidirectedEdges:
    """Test bidirected edge generation from common ancestors."""

    def test_common_ancestor_same_distance(self):
        """H -> 1, H -> 2, hide {H} => 1 <-> 2 with diff=0"""
        graph = {
            3: {1: {1: {1}}, 2: {1: {1}}},
        }
        hidden = {3}
        expected = {1: {2: {2: {0}}}}

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, expected)

    def test_common_ancestor_different_distance(self):
        """H -> A -> 1, H -> 2, hide {H, A} => 1 <-> 2 with diff=1"""
        graph = {
            3: {4: {1: {1}}, 2: {1: {1}}},
            4: {1: {1: {1}}},
        }
        hidden = {3, 4}
        expected = {1: {2: {2: {1}}}}

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, expected)

    def test_self_bidirected(self):
        """H -> A -> X, H -> X, hide {H, A} => X <-> X with diff=1"""
        graph = {
            'H': {'A': {1: {1}}, 'X': {1: {1}}},
            'A': {'X': {1: {1}}},
        }
        hidden = {'H', 'A'}
        expected = {'X': {'X': {2: {1}}}}

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, expected)


class TestLoops:
    """Test hidden nodes with self-loops."""

    def test_hidden_self_loop(self):
        """1 -> 2 -> 3 with 2 -> 2, hide {2} => 1 -> 3 with lags {2,3,4,...}"""
        graph = {
            1: {2: {1: {1}}},
            2: {2: {1: {1}}, 3: {1: {1}}},
        }
        hidden = {2}
        maxlag = 10

        expected_dir_lags = set(range(2, maxlag + 1))
        expected_bidir_lags = set(range(1, maxlag + 1))
        expected = {
            1: {3: {1: expected_dir_lags}},
            3: {3: {2: expected_bidir_lags}}
        }

        result = hide_nodes_z3(graph, hidden, maxlag=maxlag)
        assert graphs_equal(result, expected)


class TestRings:
    """Test ring/cycle topologies."""

    def test_small_ring(self):
        """Ring 1 -> 2 -> 3 -> 4 -> 5 -> 1, hide {3}"""
        graph = {
            1: {2: {1: {1}}},
            2: {3: {1: {1}}},
            3: {4: {1: {1}}},
            4: {5: {1: {1}}},
            5: {1: {1: {1}}},
        }
        hidden = {3}
        expected = {
            1: {2: {1: {1}}},
            2: {4: {1: {2}}},
            4: {5: {1: {1}}},
            5: {1: {1: {1}}},
        }

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, expected)


class TestUserGraph:
    """Test the specific graph from user request."""

    def test_user_graph_with_self_loop(self):
        """User's graph: 5 -> {6,7}, 6 <-> 7 cycle, 7 -> 8, 5 -> 5"""
        graph = {
            5: {6: {1: {1}}, 7: {1: {1}}, 5: {1: {1}}},
            6: {7: {1: {1}}},
            7: {8: {1: {1}}, 6: {1: {1}}},
            8: {}
        }
        hidden = {6, 7}
        maxlag = 17

        result = hide_nodes_z3(graph, hidden, maxlag=maxlag)

        # Check directed 5 -> 5 preserved
        assert 5 in result
        assert 5 in result[5]
        assert 1 in result[5][5]
        assert result[5][5][1] == {1}

        # Check directed 5 -> 8 has consecutive lags {2, 3, ..., 17}
        assert 8 in result[5]
        assert 1 in result[5][8]
        assert result[5][8][1] == set(range(2, maxlag + 1))

        # Check bidirected self-loop 8 <-> 8 has even diffs {2, 4, ..., 16}
        assert 8 in result
        assert 8 in result[8]
        assert 2 in result[8][8]
        assert result[8][8][2] == {2, 4, 6, 8, 10, 12, 14, 16}

    def test_user_graph_without_self_loop(self):
        """User's graph without the 5 -> 5 self-loop"""
        graph = {
            5: {6: {1: {1}}, 7: {1: {1}}},
            6: {7: {1: {1}}},
            7: {8: {1: {1}}, 6: {1: {1}}},
            8: {}
        }
        hidden = {6, 7}
        maxlag = 17

        result = hide_nodes_z3(graph, hidden, maxlag=maxlag)

        # Check directed 5 -> 8 has consecutive lags
        assert 5 in result
        assert 8 in result[5]
        assert 1 in result[5][8]
        assert result[5][8][1] == set(range(2, maxlag + 1))

        # Check bidirected self-loop 8 <-> 8
        assert 8 in result
        assert 8 in result[8]
        assert 2 in result[8][8]
        assert result[8][8][2] == {2, 4, 6, 8, 10, 12, 14, 16}


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_no_hidden_nodes(self):
        """Graph with no hidden nodes should be unchanged."""
        graph = {
            1: {2: {1: {1}}},
            2: {3: {1: {1}}},
        }
        hidden = set()

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert graphs_equal(result, graph)

    def test_all_hidden(self):
        """All nodes hidden => empty result."""
        graph = {
            1: {2: {1: {1}}},
        }
        hidden = {1, 2}

        result = hide_nodes_z3(graph, hidden, maxlag=10)
        assert result == {} or all(not nbrs for nbrs in result.values())

    def test_maxlag_truncation(self):
        """Lags should be truncated at maxlag."""
        graph = {
            1: {2: {1: {1}}},
            2: {2: {1: {1}}, 3: {1: {1}}},
        }
        hidden = {2}
        maxlag = 5

        result = hide_nodes_z3(graph, hidden, maxlag=maxlag)

        # Should only have lags up to maxlag
        assert 1 in result
        for v, ed in result.get(1, {}).items():
            for et, lags in ed.items():
                assert all(lag <= maxlag for lag in lags)

        # Specifically, 1->3 should have {2,3,4,5} not beyond
        assert 3 in result[1]
        assert result[1][3][1] == {2, 3, 4, 5}
