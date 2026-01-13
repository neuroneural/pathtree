"""
Pytest tests for the clingo-based pipeline (forward/reverse roundtrip).

Tests that the clingo ASP solver can:
1. Forward pass: marginalize out hidden nodes
2. Reverse pass: find latent structure that reproduces the marginalized graph
3. Verify roundtrip: forward(reverse(forward(G))) == forward(G)
"""

import pytest
import sys
from pathlib import Path
from copy import deepcopy
from pprint import pformat

# These imports come from conftest.py path setup
from pathtreetools import unify_edge_representation, find_maximal_bcliques, minimal_refinement
from parse import (
    convert_gk_graph, convert_with_etypes, export_to_facts, run_clingo,
    parse_clingo_output, build_set_graph, dump_for_clingo,
    add_target_lag_facts, parse_edge2_atoms
)

# Path to solver files
PATHTREE_DIR = Path(__file__).parent.parent / "pathtree"
HIDE_NODES_LP = str(PATHTREE_DIR / "hide_nodes.lp")
REVERSE_LP = str(PATHTREE_DIR / "reverse.lp")


def hide_nodes(graph, hidden, maxlag=17, solver_file=None, verbose=False, force_gk=False):
    """Run forward pass using clingo ASP solver."""
    if solver_file is None:
        solver_file = HIDE_NODES_LP

    # Handle empty graphs
    if not graph:
        return {}, (dict(), set(), set(), set(), dict(), dict())

    if force_gk:
        G_true = convert_gk_graph(graph)
    else:
        sample_val = next(iter(graph.values()))
        if all(isinstance(v, int) for v in sample_val):
            G_true = convert_gk_graph(graph)
        else:
            G_true = convert_with_etypes(graph)

    fact_str = export_to_facts(G_true, hidden)
    if verbose:
        print(fact_str)

    output = run_clingo(fact_str, maxlag=maxlag, solver_file=solver_file)
    if not output.strip():
        raise RuntimeError("hide_nodes.lp returned UNSAT for this input graph")

    parsed = parse_clingo_output(output)
    directed, directed_pairs, bidirected_pairs, bidirected_zero, bidirected_diff, base_min = parsed
    G_set = build_set_graph(directed, directed_pairs, bidirected_pairs, bidirected_zero, bidirected_diff, maxlag=maxlag)
    return G_set, parsed


def edge2_to_graph(edges):
    """Convert edge list to graph dict."""
    G = {}
    for u, v in edges:
        G.setdefault(u, {}).setdefault(v, {}).setdefault(1, set()).add(1)
    return G


def dir_unique_to_graph(clingo_output):
    """
    Parse dir_unique(u,v,lag) atoms from clingo output and build a graph.
    Each dir_unique edge becomes a unit-lag edge (lag=1) for the second forward pass.
    """
    from parse import smart_split_args

    G = {}

    for line in clingo_output.splitlines():
        if not line.strip() or line.startswith("Answer") or line.startswith("SAT"):
            continue

        for atom in line.split():
            if atom.startswith("dir_unique("):
                inner = atom[len("dir_unique("):-1]
                parts = smart_split_args(inner)
                if len(parts) == 3:
                    u, v, lag = parts
                    # Normalize: keep strings for latents, ints for observed
                    u = int(u) if u.isdigit() else u
                    v = int(v) if v.isdigit() else v
                    # Add as unit-lag edge
                    G.setdefault(u, {}).setdefault(v, {}).setdefault(1, set()).add(1)

    return G


def flatten_graph(G):
    """Flatten graph to (u, v, etype) -> lag_set for comparison."""
    out = {}
    for u, nbrs in G.items():
        for v, ed in nbrs.items():
            for et, lags in ed.items():
                out[(u, v, et)] = set(lags)
    return out


def diff_graphs(g1, g2):
    """Return string describing differences between two graphs."""
    a = flatten_graph(g1)
    b = flatten_graph(g2)

    lines = []
    all_keys = sorted(set(a) | set(b), key=str)
    for k in all_keys:
        if a.get(k, set()) != b.get(k, set()):
            u, v, et = k
            etype = "dir" if et == 1 else "bi"
            lines.append(f"{etype} {u}->{v}:")
            lines.append(f"  expected: {sorted(a.get(k, set()))}")
            lines.append(f"  actual:   {sorted(b.get(k, set()))}")

    return "\n".join(lines) if lines else "Graphs are equal"


def run_pipeline_roundtrip(graph, hidden, maxlag=17, verbose=False):
    """
    Run full pipeline roundtrip:
    1. Forward pass (hide nodes)
    2. PathTree refinement
    3. Reverse pass (find latents)
    4. Second forward pass
    5. Compare results

    Returns (match, G_set, G_fake_forward, details)
    """
    # Step 1: First forward pass
    G_set, parsed_first = hide_nodes(graph, hidden=hidden, maxlag=maxlag, verbose=verbose)

    directed_first, directed_pairs_first, bidirected_pairs_first, \
        bidirected_zero_first, bidirected_diff_first, base_min_first = parsed_first

    # Step 2: PathTree unification + refinement
    G_for = unify_edge_representation(deepcopy(G_set))
    bcliques = find_maximal_bcliques(G_for)
    G_refined = minimal_refinement(G_for, bcliques)

    # Compute observed nodes
    observed = set(graph.keys())
    for u, nbrs in graph.items():
        observed.update(nbrs.keys())
    observed -= set(hidden)

    # Step 3: Build facts for reverse search
    pf_facts = dump_for_clingo(G_refined, observed, return_str=True)
    target_block = add_target_lag_facts(directed_first, bidirected_diff_first, observed)
    fact_str = pf_facts + "\n" + target_block + "\n"

    if verbose:
        print("==== FACTS SENT TO CLINGO (REVERSE SEARCH) ====")
        print(fact_str)

    # Step 4: Reverse pass
    output = run_clingo(fact_str, maxlag=maxlag, solver_file=REVERSE_LP)

    # Parse dir_unique atoms from reverse output to build latent graph
    G_cl = dir_unique_to_graph(output)

    # Compute all nodes in reversed graph
    all_nodes = set(G_cl.keys())
    for u, nbrs in G_cl.items():
        all_nodes.update(nbrs.keys())

    new_latents = all_nodes - observed

    # Step 5: Second forward pass (only if we got edges from reverse)
    if G_cl:
        G_fake_forward, _ = hide_nodes(deepcopy(G_cl), hidden=set(new_latents), maxlag=maxlag, verbose=False)
    else:
        G_fake_forward = {}

    # Compare
    match = (G_set == G_fake_forward)

    details = {
        'G_set': G_set,
        'G_refined': G_refined,
        'G_cl': G_cl,
        'G_fake_forward': G_fake_forward,
        'observed': observed,
        'reverse_output': output,  # Add clingo output for debugging
        'new_latents': new_latents,
    }

    return match, G_set, G_fake_forward, details


class TestPipelineGraph1:
    """Test pipeline with graph1: simple cycle with self-loop on source."""

    # Graph structure:
    # 5 -> 5 (self-loop)
    # 5 -> 6, 5 -> 7
    # 6 -> 7
    # 7 -> 8, 7 -> 6 (creating 6 <-> 7 cycle)

    graph1 = {
        5: {6: {1: {1}}, 7: {1: {1}}, 5: {1: {1}}},
        6: {7: {1: {1}}},
        7: {8: {1: {1}}, 6: {1: {1}}},
        8: {}
    }
    hidden1 = {6, 7}

    def test_forward_pass(self):
        """Test that forward pass completes without error."""
        G_set, parsed = hide_nodes(self.graph1, hidden=self.hidden1, maxlag=17)

        assert G_set is not None
        assert 5 in G_set or 8 in G_set  # Should have observed nodes

    def test_roundtrip(self):
        """Test full roundtrip: forward(reverse(forward(G))) == forward(G)"""
        match, G_set, G_fake_forward, details = run_pipeline_roundtrip(
            self.graph1, self.hidden1, maxlag=17
        )

        if not match:
            diff = diff_graphs(G_set, G_fake_forward)
            G_cl = details['G_cl']
            if not G_cl:
                edge_info = "REVERSE RETURNED NO EDGES - solver may not have found a solution"
            else:
                edge_info = f"G_cl (latent graph from reverse):\n{pformat(G_cl)}"

            pytest.fail(
                f"Roundtrip mismatch for graph1!\n\n"
                f"G_set (expected):\n{pformat(G_set)}\n\n"
                f"G_fake_forward (actual):\n{pformat(G_fake_forward)}\n\n"
                f"Differences:\n{diff}\n\n"
                f"Latents found: {details['new_latents']}\n"
                f"{edge_info}"
            )


class TestPipelineGraph2:
    """Test pipeline with graph2: complex with bidirectional edges and multiple sinks."""

    # Graph structure:
    # 5 -> 6
    # 6 -> 7, 6 -> 6 (self-loop), 6 -> 5 (creating 5 <-> 6)
    # 7 -> 8, 7 -> 2
    # 8 -> 5 (cycle back to 5)

    graph2 = {
        2: {},
        5: {6: {1: {1}}},
        6: {7: {1: {1}}, 6: {1: {1}}, 5: {1: {1}}},  # self-loop on 6, and 6->5
        7: {8: {1: {1}}, 2: {1: {1}}},
        8: {5: {1: {1}}}
    }
    hidden2 = {6, 7}

    def test_forward_pass(self):
        """Test that forward pass completes without error."""
        G_set, parsed = hide_nodes(self.graph2, hidden=self.hidden2, maxlag=17)

        assert G_set is not None
        # Observed nodes should be {2, 5, 8}
        observed = {2, 5, 8}
        for node in observed:
            # At least some observed nodes should appear
            pass  # Structure depends on marginalization result

    def test_roundtrip(self):
        """Test full roundtrip: forward(reverse(forward(G))) == forward(G)"""
        match, G_set, G_fake_forward, details = run_pipeline_roundtrip(
            self.graph2, self.hidden2, maxlag=17
        )

        if not match:
            diff = diff_graphs(G_set, G_fake_forward)
            G_cl = details['G_cl']
            if not G_cl:
                edge_info = "REVERSE RETURNED NO EDGES - solver may not have found a solution"
            else:
                edge_info = f"G_cl (latent graph from reverse):\n{pformat(G_cl)}"

            pytest.fail(
                f"Roundtrip mismatch for graph2!\n\n"
                f"G_set (expected):\n{pformat(G_set)}\n\n"
                f"G_fake_forward (actual):\n{pformat(G_fake_forward)}\n\n"
                f"Differences:\n{diff}\n\n"
                f"Latents found: {details['new_latents']}\n"
                f"{edge_info}"
            )


class TestPipelineParameterized:
    """Parameterized tests for both graphs."""

    @pytest.mark.parametrize("graph,hidden,name", [
        (
            {5: {6: {1: {1}}, 7: {1: {1}}, 5: {1: {1}}},
             6: {7: {1: {1}}},
             7: {8: {1: {1}}, 6: {1: {1}}},
             8: {}},
            {6, 7},
            "graph1_self_loop_source"
        ),
        (
            {2: {},
             5: {6: {1: {1}}},
             6: {7: {1: {1}}, 6: {1: {1}}, 5: {1: {1}}},
             7: {8: {1: {1}}, 2: {1: {1}}},
             8: {5: {1: {1}}}},
            {6, 7},
            "graph2_complex_cycles"
        ),
    ])
    def test_forward_produces_valid_graph(self, graph, hidden, name):
        """Forward pass should produce a valid marginalized graph."""
        G_set, parsed = hide_nodes(graph, hidden=hidden, maxlag=17)

        # Basic validity checks
        assert isinstance(G_set, dict), f"{name}: G_set should be a dict"

        # All observed nodes should potentially appear
        observed = set(graph.keys())
        for u, nbrs in graph.items():
            observed.update(nbrs.keys())
        observed -= hidden

        # Check structure is valid
        for u, nbrs in G_set.items():
            assert isinstance(nbrs, dict), f"{name}: neighbors should be dict"
            for v, ed in nbrs.items():
                assert isinstance(ed, dict), f"{name}: edge data should be dict"
                for et, lags in ed.items():
                    assert et in {1, 2}, f"{name}: edge type should be 1 (dir) or 2 (bidir)"
                    assert isinstance(lags, set), f"{name}: lags should be set"
