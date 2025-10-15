import subprocess
from pathtree import PathTree

def convert_gk_graph(G):
    """GraphKit’s str→list repr → int→int→{1:{1}} format.
    Preserve symbolic latent IDs like h0(...), p(...).
    """
    def norm(x):
        if isinstance(x, int): 
            return x
        if isinstance(x, str) and x.isdigit():
            return int(x)
        return str(x)  # keep symbolic as string

    return {
        norm(u): {norm(v): {1: {1}} for v in G[u]}
        for u in G
    }

def convert_with_etypes(G):
    """
    Normalize a graph that's already in etype form (or close).
    Ensures every edge maps to {etype: set[int]}.
    Accepted payloads:
      - {etype: set[int]}
      - {etype: int}
      - {etype: {lag: something}} (legacy, keys must be ints)
      - plain int    (assume etype=1)
      - plain set    (assume etype=1)
    """
    out = {}
    for u, nbrs in G.items():
        out[u] = {}
        for v, ed in nbrs.items():
            out[u][v] = {}
            print(f"[convert_with_etypes] START {u}->{v}, raw={ed}")

            # Case 1: already {etype: ...}
            if isinstance(ed, dict):
                for et, lags in ed.items():
                    if isinstance(lags, set):
                        out[u][v][et] = set(lags)
                        print(f"[convert_with_etypes] normalized {u}->{v}, etype={et}, from set {lags}")

                    elif isinstance(lags, int):
                        out[u][v][et] = {lags}
                        print(f"[convert_with_etypes] normalized {u}->{v}, etype={et}, from int {lags}")

                    elif isinstance(lags, dict):
                        # Legacy style: {etype: {lag: payload}}
                        if all(isinstance(k, int) for k in lags.keys()):
                            out[u][v][et] = set(lags.keys())
                            print(f"[convert_with_etypes] normalized {u}->{v}, etype={et}, from dict keys {list(lags.keys())}")
                        else:
                            raise TypeError(
                                f"[convert_with_etypes] Unsupported dict payload for {u}->{v}, etype={et}: {lags}"
                            )

                    else:
                        raise TypeError(
                            f"[convert_with_etypes] Unexpected lag type {type(lags)} in {u}->{v}, etype={et}, val={lags}"
                        )

            # Case 2: plain int → force etype=1
            elif isinstance(ed, int):
                out[u][v][1] = {ed}
                print(f"[convert_with_etypes] int payload {ed} for {u}->{v}, forcing etype=1")

            # Case 3: plain set[int] → force etype=1
            elif isinstance(ed, set):
                if all(isinstance(x, int) for x in ed):
                    out[u][v][1] = set(ed)
                    print(f"[convert_with_etypes] set payload {ed} for {u}->{v}, forcing etype=1")
                else:
                    raise TypeError(
                        f"[convert_with_etypes] Non-int element in set payload {ed} for {u}->{v}"
                    )

            else:
                raise TypeError(f"[convert_with_etypes] Unexpected edge format {ed} for {u}->{v}")

            print(f"[convert_with_etypes] FINAL {u}->{v} -> {out[u][v]}")

    return out

def export_to_facts(G_true, hidden): 
    """Export graph into ASP facts for Clingo (hide_nodes)."""
    facts = []
    hidden = set(hidden)

    # collect ALL nodes: sources + sinks
    all_nodes = set(G_true.keys())
    for u, nbrs in G_true.items():
        all_nodes.update(nbrs.keys())

    def emit_node(x):
        if isinstance(x, int) or (isinstance(x, str) and x.isdigit()):
            return str(int(x))
        if isinstance(x, str) and "(" in x and ")" in x:
            return x  # symbolic latent like h0(5,8)
        return str(x)

    # declare nodes (deduplicated, normalized)
    seen_nodes = set()
    for v in sorted(all_nodes, key=str):
        norm_v = emit_node(v)
        if norm_v in seen_nodes:
            continue
        seen_nodes.add(norm_v)
        fact = f"node({norm_v})."
        facts.append(fact)
        print(f"[export_to_facts] DECLARE {fact}")


    # edges
    print("[export_to_facts] SNAPSHOT:", {
        u: {v: {et: sorted(list(lags)) for et, lags in ed.items()}
            for v, ed in nbrs.items()}
        for u, nbrs in G_true.items()
    })

    for u, nbrs in G_true.items():
        for v, ed in nbrs.items():
            print(f"[export_to_facts] CHECK {u}->{v} payload={ed}")

            # emit directed edges (etype=1)
            if 1 in ed:
                fact = f"edge({emit_node(u)},{emit_node(v)})."
                facts.append(fact)
                print(f"[export_to_facts] EMIT {fact} ed={ed}")

            # emit bidirected edges (etype=2)
            if 2 in ed:
                diffs = ed[2]
                # only export if both ends are observed and lag set is not just {0}
                if u not in hidden and v not in hidden and not (len(diffs) == 1 and 0 in diffs):
                    fact = f"bi_edge({emit_node(u)},{emit_node(v)})."
                    facts.append(fact)
                    print(f"[export_to_facts] EMIT {fact} ed={ed}")
                    for d in diffs:
                        if d == 0:
                            fact = f"bi_zero({emit_node(u)},{emit_node(v)})."
                            facts.append(fact)
                            print(f"[export_to_facts] EMIT {fact}")
                        else:
                            fact = f"bi_diff({emit_node(u)},{emit_node(v)},{d})."
                            facts.append(fact)
                            print(f"[export_to_facts] EMIT {fact}")
                    if diffs:
                        Dmin = min(diffs)
                        fact = f"base_min({emit_node(u)},{emit_node(v)},{Dmin})."
                        facts.append(fact)
                        print(f"[export_to_facts] EMIT {fact} (baseline distance)")
                else:
                    print(f"[export_to_facts] SKIP spurious bidirected {u}<->{v} ed={ed}")

    # observed / hidden (deduplicated & exclusive)
    all_nodes = sorted(set(all_nodes), key=str)  # ensure uniqueness again

    seen_obs = set()
    for v in all_nodes:
        norm_v = emit_node(v)
        if norm_v in seen_obs:
            continue
        if v in hidden:
            facts.append(f"hidden({norm_v}).")
        else:
            facts.append(f"observed({norm_v}).")
        seen_obs.add(norm_v)


    return "\n".join(facts)


def run_clingo(fact_str, maxlag=20, solver_path="clingo", solver_file="hide_nodes.lp"):
    with open("input.lp", "w") as f:
        f.write(fact_str)
    result = subprocess.run(
        [solver_path, solver_file, "input.lp", "--const", f"maxlag={maxlag}"],
        capture_output=True, text=True
    )

    print("STDOUT:\n", result.stdout)
    print("STDERR:\n", result.stderr)

    if result.returncode not in (0, 10, 30):
        raise RuntimeError(f"Clingo process failed with code {result.returncode}")

    if "UNSATISFIABLE" in result.stdout:
        return ""  # no model
    return result.stdout

def smart_split_args(s: str):
    """
    Split a Clingo atom argument list into top-level arguments,
    respecting parentheses (so commas inside (...) are ignored).
    Example: "h0(5,8),8,1" -> ["h0(5,8)", "8", "1"]
    """
    parts, buf, depth = [], "", 0
    for ch in s:
        if ch == ',' and depth == 0:
            parts.append(buf.strip())
            buf = ""
        else:
            if ch == '(':
                depth += 1
            elif ch == ')':
                depth -= 1
            buf += ch
    if buf:
        parts.append(buf.strip())
    return parts


def parse_clingo_output(stdout: str):
    """
    Parse clingo stdout into directed/bidirected structures.
    """
    directed = {}
    directed_pairs = set()
    bidirected_pairs = set()
    bidirected_zero = set()
    bidirected_diff = {}
    def norm(x):
        x = x.strip()
        if x.isdigit():
            return int(x)
        return x  # keep symbolic names like h0(5,8), p(...), cyc(...)

    for line in stdout.splitlines():
        if not line.strip() or line.startswith("Answer") or line.startswith("SAT") or line.startswith("Models"):
            continue

        for atom in line.split():
            if atom.startswith("latent(pu(") or atom.startswith("latent(bcyc("):
                continue  # skip helper latent chains entirely
            if atom.startswith("dir_unique("):
                inner = atom[len("dir_unique("):-1]
                parts = smart_split_args(inner)
                if len(parts) != 3:
                    raise ValueError(f"Unexpected dir_unique format: {inner}")
                u, v, L = parts
                key = (norm(u), norm(v))
                directed.setdefault(key, set()).add(int(L))

            elif atom.startswith("dir_edge("):
                inner = atom[len("dir_edge("):-1]
                parts = smart_split_args(inner)
                if len(parts) != 2:
                    raise ValueError(f"Unexpected dir_edge format: {inner}")
                u, v = parts
                directed_pairs.add((norm(u), norm(v)))

            elif atom.startswith("bi_edge("):
                inner = atom[len("bi_edge("):-1]
                u, v = smart_split_args(inner)
                bidirected_pairs.add((norm(u), norm(v)))

            elif atom.startswith("bi_zero("):
                inner = atom[len("bi_zero("):-1]
                u, v = smart_split_args(inner)
                bidirected_zero.add((norm(u), norm(v)))

            elif atom.startswith("bi_diff("):
                inner = atom[len("bi_diff("):-1]
                u, v, d = smart_split_args(inner)
                key = (norm(u), norm(v))
                bidirected_diff.setdefault(key, set()).add(int(d))

            elif atom.startswith("bi_unique("):
                inner = atom[len("bi_unique("):-1]
                parts = smart_split_args(inner)
                if len(parts) != 4:
                    raise ValueError(f"Unexpected bi_unique format: {inner}")
                # we don't actually use bi_unique payloads in set-graph building
                u, v, _, _ = parts
                bidirected_pairs.add((norm(u), norm(v)))
            # elif atom.startswith("base_min_final("):
            #     inner = atom[len("base_min_final("):-1]
            #     x, y, dmin = smart_split_args(inner)
            #     key = (norm(x), norm(y))
            #     val = int(dmin)
            #     prev = d_child.get(key)
            #     if prev is None or val < prev:
            #         d_child[key] = val




            # ignore out_node/observed/latent/hidden for graph build
            # but you could debug-log them if needed

    return directed, directed_pairs, bidirected_pairs, bidirected_zero, bidirected_diff

def dump_for_clingo(graph, observed, base_min=None, return_str=False, filename=None): 
    """
    Serialize a PathForest/lag-set graph into ASP facts for reverse.lp.
    Instead of flattening to dir_unique(U,V,L), we preserve path-tree structure:
      - root(U,V,A) for a path-tree root with preset A
      - loop(U,V,D) for each child cycle of length D
    If the edge is already a plain finite lag set, emit dir_unique(U,V,L).
    """
    facts = []

    # emit observed nodes explicitly
    for u in sorted(observed):
        facts.append(f"observed({str(u)}).")

    for u, nbrs in graph.items():
        for v, ed in nbrs.items():
            for et, raw in ed.items():

                # ===============================
                # Bidirected edges
                # ===============================
                if et == 2:
                    # canonicalize orientation so we only emit once
                    uu, vv = (u, v) if str(u) < str(v) else (v, u)
                    facts.append(f"bi_edge({uu},{vv}).")
                    if isinstance(raw, PathTree) and raw.loopset:
                        for child in raw.loopset:
                            if isinstance(child, int):
                                facts.append(f"loop({uu},{vv},{child}).")
                            elif isinstance(child, PathTree):
                                subpreset = next(iter(child.preset)) if isinstance(child.preset, set) else child.preset
                                facts.append(f"loop({uu},{vv},{subpreset}).")
                    if isinstance(raw, set):
                        for d in sorted(raw):
                            if int(d) == 0:
                                facts.append(f"bi_zero({uu},{vv}).")
                            else:
                                facts.append(f"bi_diff({uu},{vv},{int(d)}).")

                    elif isinstance(raw, list):
                        for elt in raw:
                            # plain int
                            if isinstance(elt, int):
                                if elt == 0:
                                    facts.append(f"bi_zero({uu},{vv}).")
                                else:
                                    facts.append(f"bi_diff({uu},{vv},{elt}).")

                            # PathTree with a numeric preset
                            elif isinstance(elt, PathTree):
                                preset = elt.preset
                                if isinstance(preset, set) and len(preset) == 1:
                                    preset = next(iter(preset))
                                if isinstance(preset, int):
                                    if preset == 0:
                                        facts.append(f"bi_zero({uu},{vv}).")
                                    else:
                                        facts.append(f"bi_diff({uu},{vv},{preset}).")
                                else:
                                    raise TypeError(f"Unsupported PathTree preset in bidirected: {preset}")

                            # tuple of ints
                            elif isinstance(elt, tuple) and all(isinstance(x, int) for x in elt):
                                for val in elt:
                                    if val == 0:
                                        facts.append(f"bi_zero({uu},{vv}).")
                                    else:
                                        facts.append(f"bi_diff({uu},{vv},{val}).")

                            else:
                                raise TypeError(f"Unsupported bidirected element in list: {elt}")

                    elif isinstance(raw, PathTree):
                        preset = raw.preset
                        if isinstance(preset, set) and len(preset) == 1:
                            preset = next(iter(preset))
                        if isinstance(preset, int):
                            if preset == 0:
                                facts.append(f"bi_zero({uu},{vv}).")
                            else:
                                facts.append(f"bi_diff({uu},{vv},{preset}).")
                    else:
                        raise TypeError(f"Unsupported bidirected payload: {raw}")

                    continue  # done with this edge

                # Directed edges
                if et != 1:
                    continue

                # PathTree case
                if isinstance(raw, PathTree):
                    preset = raw.preset
                    if isinstance(preset, set) and len(preset) == 1:
                        preset = next(iter(preset))
                    elif isinstance(preset, set):
                        raise ValueError(f"Unsupported multi-base preset: {preset}")

                    facts.append(f"root({u},{v},{preset}).")

                    # Emit loops (flat children only for now)
                    for child in raw.loopset:
                        if isinstance(child, int):
                            facts.append(f"loop({u},{v},{child}).")
                        elif isinstance(child, PathTree):
                            sub_str = dump_for_clingo(
                                {u: {v: {1: child}}},
                                observed,
                                return_str=True
                            )
                            facts.extend(sub_str.strip().splitlines())
                        else:
                            raise TypeError(f"Unsupported loopset element: {child}")

                # finite lag set
                elif isinstance(raw, set):
                    for L in sorted(raw):
                        facts.append(f"dir_unique({u},{v},{L}).")

                # list of PathTrees or ints
                elif isinstance(raw, list):
                    for elt in raw:
                        if isinstance(elt, PathTree):
                            sub_str = dump_for_clingo(
                                {u: {v: {1: elt}}},
                                observed,
                                return_str=True
                            )
                            facts.extend(sub_str.strip().splitlines())
                        elif isinstance(elt, int):
                            facts.append(f"dir_unique({u},{v},{elt}).")
                        else:
                            raise TypeError(f"Unsupported elt in list: {elt}")
                else:
                    raise TypeError(f"Unsupported edge payload: {raw}")
    if base_min:
        for (u, v), dmin in base_min.items():
            facts.append(f"base_min_final({u},{v},{dmin}).")


    out = "\n".join(facts) + "\n"

    if return_str:
        return out
    if filename:
        with open(filename, "w") as f:
            f.write(out)
    return None

def build_set_graph(directed, directed_pairs, bidirected_pairs, bidirected_zero,
                    bidirected_diff, maxlag=17):
    G = {}
    sinks = set()

    def norm(x):
        try:
            return int(x)
        except ValueError:
            return x

    # dir_unique(u,v,L)
    for (u, v), lags in directed.items():
        u, v = norm(u), norm(v)
        sinks.add(v)
        expanded = set()
        for L in lags:
            # if u == v and L > 1:
            #     expanded.update({k * L for k in range(1, maxlag // L + 1)})
            # else:
            expanded.add(L)
        G.setdefault(u, {}).setdefault(v, {})[1] = expanded
        print(f"[DEBUG build] dir_unique made directed {u}->{v} lags={sorted(expanded)}")

    # raw unit edges edge(u,v)
    for (u, v) in directed_pairs:
        u, v = norm(u), norm(v)
        sinks.add(v)
        # only add lag=1 if this pair has *no explicit dir_unique entry*
        if (u, v) not in directed:
            if u != v:
                G.setdefault(u, {}).setdefault(v, {}).setdefault(1, set()).add(1)
                print(f"[DEBUG build] edge/2 made directed {u}->{v} add lag 1")
            else:
                print(f"[DEBUG build] skipped spurious unit self-edge {u}->{v}")
        else:
            print(f"[DEBUG build] skipped redundant edge/2 for {u}->{v}")



    # bidirected edges
    for (u, v) in bidirected_pairs:
        u, v = norm(u), norm(v)
        diffs = set()
        if (u, v) in bidirected_diff:
            diffs |= bidirected_diff[(u, v)]
        if (u, v) in bidirected_zero or (v, u) in bidirected_zero:
            diffs.add(0)
        # if d_child and (u, v) in d_child:
        #     base_dist = d_child[(u, v)]
        #     diffs = {d for d in diffs if d <= maxlag - base_dist}

        print(f"[DEBUG build_set_graph] Pair {(u,v)} raw diffs={diffs}")
        if not diffs:
            print(f"[DEBUG build_set_graph] Skipping {(u,v)} (no lags)")
            continue

        G.setdefault(u, {}).setdefault(v, {})[2] = diffs
        print(f"[DEBUG build_set_graph] Added bidirected {u}<->{v} with {sorted(diffs)}")

    # final snapshot of every edge dict
    print("[DEBUG build_set_graph] Final G snapshot:")
    for u, nbrs in G.items():
        for v, ed in nbrs.items():
            print(f"  {u}->{v} : {ed}")

    print(f"[DEBUG build_set_graph] Final G keys={list(G.keys())}")
    return G
