from copy import deepcopy
from gunfolds.utils.graphkit import ringmore
from pathtreetools import unify_edge_representation, find_maximal_bcliques, minimal_refinement
from parse import *
from pprint import pprint
VERBOSE=True
    
def hide_nodes(graph, hidden={6,7,10}, maxlag=17, solver_file="hide_nodes.lp", verbose=VERBOSE, force_gk=False):
    if locals().get("force_gk"):
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
    parsed = parse_clingo_output(output)
    G_set = build_set_graph(*parsed, maxlag=maxlag)
    return G_set
def normalize_ringmore(raw):
    """Convert ringmore output {u:{v:1}} → {u:{v:{1:{1}}}}."""
    G = {}
    for u, nbrs in raw.items():
        G[u] = {}
        for v in nbrs:
            G[u][v] = {1: {1}}
    return G


def format_graph_one_line(G):
    """Compact one-line-per-node representation."""
    lines = ["{"]
    for u, nbrs in G.items():
        inner = ", ".join(f"{v}: {ed}" for v, ed in nbrs.items())
        lines.append(f"    {u}: {{{inner}}},")
    lines.append("}")
    return "\n".join(lines)


def run_ringmore_trials(num_trials=10, N=10, num_extra_edges=6, hidden={6, 7, 10}, maxlag=17):
    failures = 0
    first_fail_printed = False

    for i in range(1, num_trials + 1):
        print(f"\n========== Trial {i} ==========")
        result = ringmore_pipeline(i, N, num_extra_edges, hidden, maxlag)

        if not result["match"]:
            failures += 1
            if not first_fail_printed:
                print("\n===== FIRST FAIL CASE ORIGINAL GRAPH =====")
                print(format_graph_one_line(result["graph"]))
                print("==========================================")
                first_fail_printed = True

    print(f"\nTotal Failures: {failures}/{num_trials}")


def ringmore_pipeline(trial_idx, N=10, num_extra_edges=6, hidden={6, 7, 10}, maxlag=17):
    # Step 1: Generate and show the raw graph
    raw = ringmore(N, num_extra_edges)
    graph = normalize_ringmore(raw)

    # Step 2: First hide phase
    G_set = hide_nodes(convert_gk_graph(graph), hidden=hidden, maxlag=maxlag, verbose=False, force_gk=True)
    print("G_set (plain sets):")
    pprint(G_set)

    # Step 3: PathTree unification + refinement
    G_for = unify_edge_representation(deepcopy(G_set))
    print("G_for:")
    pprint(G_for)
    bcliques = find_maximal_bcliques(G_for)
    G_refined = minimal_refinement(G_for, bcliques)
    print("\nGraph with PathForest weights (refined):")
    pprint(G_refined)

    # Step 4: Reverse phase
    observed = set(graph.keys()) - set(hidden)
    fact_str = dump_for_clingo(G_refined, observed, return_str=True)
    print("==== FACTS SENT TO CLINGO ====")
    print(fact_str)
    output = run_clingo(fact_str, maxlag=maxlag, solver_file="reverse.lp")

    directed, directed_pairs, bidirected_pairs, bidirected_zero, bidirected_diff = parse_clingo_output(output)
    print("[DEBUG] Types:", type(directed), type(bidirected_pairs),
          type(bidirected_diff), type(bidirected_zero))
    print("[DEBUG] Lengths:", len(directed), len(bidirected_pairs),
          len(bidirected_diff), len(bidirected_zero))
    print("[DEBUG] Sample bidirected_diff:", bidirected_diff[:5] if isinstance(bidirected_diff, list) else list(bidirected_diff.items())[:5])
    G_cl = build_set_graph(directed, directed_pairs, bidirected_pairs, bidirected_zero, bidirected_diff)
    print("[DEBUG pre-normalize G_cl snapshot]:")
    for u, nbrs in G_cl.items():
        for v, ed in nbrs.items():
            print(f"  {u}->{v} : { {et: sorted(list(lags)) for et, lags in ed.items()} }")

    # Normalize keys (integers stay ints, latents stay strings)
    G_cl = {int(k) if isinstance(k, str) and k.isdigit() else k: v for k, v in G_cl.items()}
    print("[DEBUG] Raw G_cl keys before casting:", list(G_cl.keys())[:10])
    print(f"DEBUG: G_cl keys: {sorted(G_cl.keys(), key=str)}")
    print(f"DEBUG: observed nodes: {sorted(observed)}")

    # Step 5: Second hide phase
    new_latents = set(G_cl) - observed
    print(f"DEBUG: About to call hide_nodes with new_latents: {sorted(map(str, new_latents))}")
    G_fake_forward = hide_nodes(deepcopy(G_cl), list(new_latents), maxlag=maxlag, verbose=False)

    print("\nGraph after second forward pass (G_fake_forward):")
    pprint(G_fake_forward)
    print("\nCompare only lag sets:", G_set == G_fake_forward)

    match = G_set == G_fake_forward
    print(f"Trial {trial_idx}: Match = {match}")

    return {"graph": graph, "match": match}

if __name__ == "__main__":
    run_ringmore_trials(num_trials=1000, N=10, num_extra_edges=6, hidden={6, 7, 10}, maxlag=17)