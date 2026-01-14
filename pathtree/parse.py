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

            if 1 in ed:
                lags = ed[1]
                if u == v and u in hidden:
                    # hidden self-loop: always encode as loop_len, including 1
                    for d in sorted(lags):
                        facts.append(f"loop_len({emit_node(u)},{int(d)}).")
                        print(f"[export_to_facts] EMIT loop_len({emit_node(u)},{int(d)}).")
                    continue  # IMPORTANT: do not also emit edge(u,u)
                else:
                    facts.append(f"edge({emit_node(u)},{emit_node(v)}).")
                    print(f"[export_to_facts] EMIT edge({emit_node(u)},{emit_node(v)}). ed={ed}")


            # emit bidirected edges (etype=2)
            if 2 in ed:
                diffs = ed[2]
                # export all observed–observed bidirected edges unless it is a pure self-loop {0}
                if u not in hidden and v not in hidden and not (u == v and len(diffs) == 1 and 0 in diffs):
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
                            # NEW: mark this as an input bidirected lag
                            input_fact = f"input_bi_diff({emit_node(u)},{emit_node(v)},{d})."
                            facts.append(input_fact)
                            print(f"[export_to_facts] EMIT {input_fact} (input provenance)")
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


def run_clingo(fact_str, maxlag=20, solver_path="clingo", solver_file="hide_nodes.lp", verbose=False):
    with open("input.lp", "w") as f:
        f.write(fact_str)
    result = subprocess.run(
        [solver_path, solver_file, "input.lp", "--const", f"maxlag={maxlag}"],
        capture_output=True, text=True
    )

    if verbose:
        print("STDOUT:\n", result.stdout)
        if result.stderr.strip():
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
    base_min = {}
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
            elif atom.startswith("loop("):
                inner = atom[len("loop("):-1]
                u, v, d = smart_split_args(inner)
                nu = norm(u)
                nv = norm(v)

                # Only treat self-loops on *latent* nodes as directed loops.
                # Observed nodes are plain ints (5, 8, ...), latents are strings like "h0(5,8)", "p(5,8,2,1)".
                if nu == nv and not isinstance(nu, int):
                    key = (nu, nv)
                    directed.setdefault(key, set()).add(int(d))
                    print(f"[DEBUG parse] captured latent self-loop {u}->{v} lag={d}")
                else:
                    # For observed nodes (nu is int) we ignore loop/3 when building the set-graph.
                    # The loop info is only for reverse's internal PathTree structure, not a real dir self-edge.
                    print(f"[DEBUG parse] ignoring observed self-loop loop({u},{v},{d}) for directed graph")

            elif atom.startswith("base_min_final("):
                inner = atom[len("base_min_final("):-1]
                u, v, d = smart_split_args(inner)
                key = (norm(u), norm(v))
                base_min[key] = int(d)

    return directed, directed_pairs, bidirected_pairs, bidirected_zero, bidirected_diff, base_min

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

                if et == 2:
                    uu, vv = u, v

                    # Collect all lags in this bidirected edge (whatever its type)
                    lags = set()
                    if isinstance(raw, set):
                        lags |= raw
                    elif isinstance(raw, PathTree):
                        if isinstance(raw.preset, set):
                            lags |= raw.preset
                        elif isinstance(raw.preset, int):
                            lags.add(raw.preset)
                    elif isinstance(raw, list):
                        for elt in raw:
                            if isinstance(elt, int):
                                lags.add(elt)
                            elif isinstance(elt, tuple):
                                lags.update(elt)
                            elif isinstance(elt, PathTree):
                                if isinstance(elt.preset, set):
                                    lags |= elt.preset
                                elif isinstance(elt.preset, int):
                                    lags.add(elt.preset)
                    # emit loop() facts for bidirected self-loop PathTrees ---
                    if isinstance(raw, PathTree):
                        for loop in sorted(raw.loopset):
                            if isinstance(loop, int):
                                facts.append(f"loop({uu},{vv},{loop}).")
                    elif isinstance(raw, list):
                        for elt in raw:
                            if isinstance(elt, PathTree):
                                for loop in sorted(elt.loopset):
                                    if isinstance(loop, int):
                                        facts.append(f"loop({uu},{vv},{loop}).")

                    has_zero = 0 in lags

                    # Emit bi_edge facts:
                    # - Always emit u→v
                    # - Only emit v→u if zero exists in lag set (instantaneous)
                    facts.append(f"bi_edge({uu},{vv}).")
                    if uu != vv and has_zero:
                        facts.append(f"bi_edge({vv},{uu}).")

                    # Emit per-lag facts:
                    # - Mirror only zeros
                    for d in sorted(lags):
                        if d == 0:
                            facts.append(f"bi_zero({uu},{vv}).")
                            if uu != vv:
                                facts.append(f"bi_zero({vv},{uu}).")
                        else:
                            facts.append(f"bi_diff({uu},{vv},{d}).")

                    continue  # done with this edge


                # Directed edges (PathTree / PathForest)
                if et != 1:
                    continue

                # Case 1: PathForest (list of PathTrees)
                if isinstance(raw, list):
                    # Group PathTrees that share the same parent edge (u,v)
                    print(f"[DEBUG PathForest check] Edge ({u}->{v}) has {len(raw)} PathTrees:")
                    for pt in raw:
                        if isinstance(pt, PathTree):
                            # Extract numeric preset
                            print(f"   preset={pt.preset}, loopset={pt.loopset}")
                            if isinstance(pt.preset, set):
                                if len(pt.preset) != 1:
                                    raise ValueError(f"Multi-base preset not supported: {pt.preset}")
                                preset_val = next(iter(pt.preset))
                            else:
                                preset_val = pt.preset

                            # Emit root fact
                            facts.append(f"root({u},{v},{int(preset_val)}).")

                            # Emit each loop fact from loopset
                            for loop in sorted(pt.loopset):
                                if isinstance(loop, int):
                                    facts.append(f"loop({u},{v},{loop}).")
                                elif isinstance(loop, PathTree):
                                    sub = dump_for_clingo({u: {v: {1: loop}}}, observed, return_str=True)
                                    facts.extend(sub.strip().splitlines())
                                else:
                                    raise TypeError(f"Unsupported loopset element {loop!r}")
                        elif isinstance(pt, int):
                            # numeric lag
                            facts.append(f"dir_unique({u},{v},{pt}).")
                        else:
                            raise TypeError(f"Unsupported element in list for directed edge ({u},{v}): {type(pt)}")

                # Case 2: Single PathTree
                elif isinstance(raw, PathTree):
                    if isinstance(raw.preset, set):
                        if len(raw.preset) != 1:
                            raise ValueError(f"Multi-base preset not supported: {raw.preset}")
                        preset_val = next(iter(raw.preset))
                    else:
                        preset_val = raw.preset

                    # special case: preset=1 with loopset={1} encodes dual path
                    if preset_val == 1 and set(raw.loopset) == {1}:
                        # normal direct edge
                        facts.append(f"root({u},{v},1).")
                        # synthetic delayed path
                        facts.append(f"root({u},{v},2).")
                        # loop belongs to the delayed branch (A=2), not A=1
                        facts.append(f"loop({u},{v},1).")
                        continue

                    # --- normal single-root case ---
                    facts.append(f"root({u},{v},{int(preset_val)}).")


                    for loop in sorted(raw.loopset):
                        if isinstance(loop, int):
                            facts.append(f"loop({u},{v},{loop}).")
                        elif isinstance(loop, PathTree):
                            sub = dump_for_clingo({u: {v: {1: loop}}}, observed, return_str=True)
                            facts.extend(sub.strip().splitlines())
                        else:
                            raise TypeError(f"Unsupported loopset element {loop!r}")

                # Case 3: Finite lag set
                elif isinstance(raw, set):
                    for L in sorted(raw):
                        facts.append(f"dir_unique({u},{v},{L}).")

                else:
                    raise TypeError(f"Unsupported directed payload for edge ({u},{v}): {type(raw)}")

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
    def norm(x):
        try:
            return int(x)
        except ValueError:
            return x

    # dir_unique(u,v,L)
    for (u, v), lags in directed.items():
        u, v = norm(u), norm(v)

        # skip redundant directed self-loops that already belong to a bidirected pair
        if u == v and (u, v) in bidirected_pairs:
            print(f"[DEBUG build] keep self-loop dir_unique {u}->{v} even though bidirected exists")
        expanded = set()
        for L in lags:
            expanded.add(L)
        G.setdefault(u, {}).setdefault(v, {})[1] = expanded
        print(f"[DEBUG build] dir_unique made directed {u}->{v} lags={sorted(expanded)}")


    # raw unit edges edge(u,v)
    for (u, v) in directed_pairs:
        u, v = norm(u), norm(v)
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

        print(f"[DEBUG build_set_graph] Pair {(u,v)} raw diffs={diffs}")
        if not diffs:
            print(f"[DEBUG build_set_graph] Skipping {(u,v)} (no lags)")
            continue
        # keep directed and bidirected self-loops separate
        if u == v and 0 in diffs and len(diffs) == 1:
            print(f"[DEBUG build_set_graph] Pure zero-lag self-loop {u}<->{v} — skip bidir aggregation")
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

def add_target_lag_facts(directed, bidirected_diff, observed):
    lines = []
    # directed lagsets from first hide_nodes run
    for (u, v), lags in directed.items():
        if u in observed and v in observed:
            for L in sorted(lags):
                lines.append(f"target_dir({u},{v},{L}).")

    # bidirected lagsets
    for (u, v), diffs in bidirected_diff.items():
        if u in observed and v in observed:
            for D in sorted(diffs):
                lines.append(f"target_bi({u},{v},{D}).")

    return "\n".join(lines)

def _split_top_level_comma(s: str) -> tuple[str, str]:
    depth = 0
    for i, ch in enumerate(s):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 0:
            left = s[:i].strip()
            right = s[i+1:].strip()
            return left, right
    raise ValueError(f"No top-level comma found in: {s!r}")


def parse_edge2_atoms(clingo_out: str):
    """
    Parse edge/2 atoms from clingo output, supporting function terms like aux(1).
    Returns list of (u, v) where u,v are ints if purely numeric, else strings.
    """
    edges = []

    # Find the Answer line(s) and parse tokens from them.
    # This is simple and works with your current outputs.
    for line in clingo_out.splitlines():
        line = line.strip()
        if not line or line.startswith("Answer:") or line.startswith("Optimization:"):
            continue
        # Heuristic: lines containing atoms have "edge(" tokens
        if "edge(" not in line:
            continue

        tokens = line.split()
        for tok in tokens:
            if not tok.startswith("edge("):
                continue
            if not tok.endswith(")"):
                # sometimes clingo prints atoms without trailing ".", but should end with ")"
                # if you have "edge(...)." then strip the final "."
                tok = tok.rstrip(".")
            inner = tok[len("edge("):-1]  # remove "edge(" and trailing ")"
            u_raw, v_raw = _split_top_level_comma(inner)

            def cast(x: str):
                return int(x) if x.isdigit() else x

            edges.append((cast(u_raw), cast(v_raw)))

    return edges
