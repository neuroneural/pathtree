import sys
import heapq
sys.path.append('./tools/')
from copy import deepcopy
from itertools import combinations
from math import gcd
from matplotlib.cbook import flatten

from ortools.constraint_solver import pywrapcp
from pathtree import PathTree
import numpy as np
from sortedcontainers import SortedDict

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

def forest_to_set(forest, cap=5):
    """
    Collapse a PathForest into the first `cap` integer lengths,
    allowing loops to repeat indefinitely but stopping once cap are found.
    """
    # 1) Normalize into PathTree roots
    roots = to_path_forest(forest)
    # split any multi-base root into one root per base
    expanded = []
    for r in roots:
        if isinstance(r.preset, set) and len(r.preset) > 1:
            for b in r.preset:
                # shallow-copy the loopset so every clone shares it
                expanded.append(PathTree(preset={b}, loopset=r.loopset.copy()))
        else:
            expanded.append(r)
    roots = expanded
    # 2) Min-heap of (total_length, node)
    heap = []
    for r in roots:
        base = next(iter(r.preset)) if isinstance(r.preset, set) else r.preset
        heapq.heappush(heap, (base, r))

    result = set()

    # 3) Pop smallest total, record it, then re-enqueue loops on THAT node
    while heap and (cap is None or len(result) < cap):
        total, node = heapq.heappop(heap)
        if total in result:
            continue
        result.add(total)

        # re-enqueue loops on the same node
        for loop in node.loopset:
            if isinstance(loop, PathTree):
                # a nested loop‐tree
                L = next(iter(loop.preset)) if isinstance(loop.preset, set) else loop.preset
                heapq.heappush(heap, (total + L, node))
            elif isinstance(loop, int):
                # a flat self‐loop of length `loop`
                heapq.heappush(heap, (total + loop, node))
            else:
                # just in case someone stuffed something else in there
                try:
                    # if it's a singleton set, etc.
                    L = int(loop)
                    heapq.heappush(heap, (total + L, node))
                except:
                    raise TypeError(f"Unsupported loop item {loop!r}")

    return result

def extends(pt, lag: int, Eij: set[int]) -> bool:
    new_pt = deepcopy(pt)
    new_pt.add_loop(lag)
    # collect the first few delays (default cap=5 is enough to see any extras)
    all_lags = forest_to_set(new_pt)  
    # we only allow this loop if it doesn't introduce _any_ new delays
    return all_lags.issubset(Eij)

def refine_edge(pt: PathTree, Eij: set[int]) -> PathTree:
    """
    Starting from pt.preset = base, add exactly the deltas needed so that
    every e in Eij can be generated (i.e. e = base + sum of loops).
    """
    # 1) extract base
    if isinstance(pt.preset, set):
        assert len(pt.preset) == 1
        base = next(iter(pt.preset))
    else:
        base = pt.preset

    # 2) for each observed lag > base, try exactly delta = e-base
    for e in sorted(Eij):
        if e <= base:
            continue
        delta = e - base
        if extends(pt, delta, Eij):
            pt.add_loop(delta)
    generated = forest_to_set(pt)
    if generated != Eij:
        # fallback: a “flat” forest that exactly matches Eij, no loops
        return PathTree(preset=Eij, loopset=set())

    return pt

def extract_bidirected(G):
    next_latent = max(G) + 1
    seen = set()
    for u in list(G):
        for v, edges in list(G[u].items()):
            if 2 in edges and (u, v) not in seen:
                seen.add((u,v)); seen.add((v,u))
                for lag in edges[2]:
                    H = next_latent; next_latent += 1
                    G.setdefault(u, {})[H] = {1: {lag}}
                    G.setdefault(H, {})[v] = {1: {lag}}
                # remove both sides of the bi-directed link
                del G[u][v][2]
                if not G[u][v]: del G[u][v]
                if 2 in G.get(v, {}).get(u, {}):
                    del G[v][u][2]
                    if not G[v][u]: del G[v][u]
    return G

def forward(graph, observed_nodes, original_nodes):
    from latents import hide_node

    G = deepcopy(graph)
    # split latents so injected come first
    all_latents      = [n for n in G if n not in observed_nodes]
    injected_latents = [n for n in all_latents if n not in original_nodes]
    original_latents = [n for n in all_latents if n in     original_nodes]

    for H in injected_latents + original_latents:
        # injected ⇒ zero_incoming=True, original ⇒ zero_incoming=False
        zero_in = (H not in original_nodes)
        G = hide_node(G, H, zero_incoming=zero_in)

        # now recompress every edge exactly as before
        for u in list(G):
            for v in list(G[u]):
                forest   = G[u][v].get(1) or G[u][v].get(2)
                new_lags = forest_to_set(forest)
                et       = 1 if 1 in G[u][v] else 2
                G[u][v]  = {et: new_lags}

    return G

def full_forward(raw_graph, observed_nodes):
    # 1) record which nodes are “real”
    original = set(raw_graph.keys())

    # 2) peel off bidirected edges
    G = deepcopy(raw_graph)
    G = extract_bidirected(G)

    # 3) wrap whatever remains into PathTree forests (or keep existing ones)
    for u in list(G):
        for v in list(G[u]):
            et      = next(iter(G[u][v]))
            dat     = G[u][v][et]
            if isinstance(dat, PathTree):
                forest = [dat]
            elif isinstance(dat, list):
                forest = dat
            else:
                forest = to_path_forest(dat)
            G[u][v] = {et: forest}

    # 4) **one** backward pass is enough to remove *all* latents
    G = forward(G, observed_nodes, original)

    # 5) make sure observed nodes still appear
    for o in observed_nodes:
        G.setdefault(o, {})

    for u in list(G):
        for v in list(G[u]):
            val = G[u][v].get(1) or G[u][v].get(2)
            if isinstance(val, set):
                lags = val
            else:
                forest = val if isinstance(val, list) else [val]
                lags = forest_to_set(forest)
            G[u][v] = {1: lags}
    return G

from math import gcd
from copy import deepcopy
from pathtreetools import forest_to_set
from pathtree import PathTree          # only needed for type check

def decompress_to_unit_graph(obs_graph):
    """
    1) Collapse any PathTree to an integer lag-set.
    2) Pull each bidirected edge u<->v into a latent H.
    3) Detect arithmetic-progression directed edges and encode them
       with a single latent + self-loop (kept unexpanded).
    4) Expand every remaining integer-lag edge to unit-lag chains.
    5) Ensure all observed nodes still appear.
    """

    # ── Step 0: normalise every edge to {etype: set(int)} ──────────────
    G0 = {}
    for u, nbrs in obs_graph.items():
        G0[u] = {}
        for v, ed in nbrs.items():
            collapsed = {}
            for etype, raw in ed.items():
                if isinstance(raw, set):
                    lags = raw
                else:                           # list[PathTree] | PathTree
                    forest = raw if isinstance(raw, list) else [raw]
                    lags   = forest_to_set(forest)
                collapsed[etype] = set(lags)
            G0[u][v] = collapsed

    # ── Step 1: move bidirected edges → new latents ───────────────────
    G1 = {}
    next_latent   = max(G0) + 1
    bidir_latents = set()

    seen_bidir = set()          # NEW  unordered pairs already handled

    for u, nbrs in G0.items():
        for v, ed in nbrs.items():

            # (A) keep directed edge
            if 1 in ed:
                G1.setdefault(u, {})[v] = {1: set(ed[1])}

            # (B) explode each bidirected lag
            if 2 in ed:
                pair = tuple(sorted((u, v)))   # {u,v} as an unordered key
                if pair in seen_bidir:
                    continue                   # already processed this pair
                seen_bidir.add(pair)

                for L in sorted(ed[2]):        # keep *all* bidirected lags
                    H = next_latent; next_latent += 1
                    bidir_latents.add(H)

                    # H → u   (lag 1)
                    G1.setdefault(H, {})[u] = {1: {1}}
                    # H → v   (lag L+1)   (integer lag, no chain yet)
                    G1[H][v] = {1: {L + 1}}

    # ── Step 1½: compress A-P lag-sets into single latent + self-loop ──
    G2 = deepcopy(G1)        # we’ll delete / add edges inside the loop
    for u, nbrs in list(G1.items()):
        for v, ed in list(nbrs.items()):
            if u == v or 1 not in ed:
                continue
            L = ed[1]
            if len(L) < 2:
                continue
            a = min(L)
            d = gcd(*[x - a for x in L])
            if d == 0 or L != {a + k*d for k in range(len(L))}:
                continue        # not an exact A-P → leave unchanged
            if a == 1:
                continue
            # ----- rewrite this edge -----------------------------------
            del G2[u][v]
            if not G2[u]:            # clean empties
                del G2[u]
            # (1) u  ->  H   at lag 1
            H = next_latent; next_latent += 1
            # H -> u   (lag 1)
            G2.setdefault(u, {})[H] = {1: {1}}
            # (2) self-loop on H  (lag = d)
            G2.setdefault(H, {})[H] = {1: {d}}
            # H -> v   (lag a)  (expand to chain of (a-1) unit lags)
            prev = H
            if a == 1:
                # special case: min lag was 1 (shouldn’t actually happen here)
                G2.setdefault(prev, {})[v] = {1: {1}}
            else:
                # need (a-1) unit edges *including* the final one to v
                for _ in range(a - 2):          # <-- a-2 intermediates
                    K = next_latent; next_latent += 1
                    G2.setdefault(prev, {})[K] = {1: {1}}
                    prev = K
                G2.setdefault(prev, {})[v] = {1: {1}}   # final edge

    # directed-only graph to expand
    Gdir = G2

    # ── Step 2: expand *non-self* integer lags to unit chains ──────────
    G_unit = {}
    for u, nbrs in Gdir.items():
        G_unit.setdefault(u, {})
                # ▸ 1.  Skip the special latents that came from bidirected edges
        if u in bidir_latents:
            # just copy their integer-lag edges verbatim
            for v, ed in nbrs.items():
                G_unit[u][v] = {1: set(ed[1])}
            continue          # go straight to the next u
        # NEW: remember the helper nodes we may create for each lag length
        # **only for this parent u**
        chain_cache = {}          #  key = L   value = list of helper nodes

        for v, ed in nbrs.items():
            lags = ed[1]

            for L in sorted(lags):

                # (a) keep self-loops as they are
                if u == v:
                    G_unit[u]\
                        .setdefault(u, {})\
                        .setdefault(1, set())\
                        .add(L)
                    continue

                # (b) unit lag – nothing to expand
                if L == 1:
                    G_unit[u].setdefault(v, {})[1] = {1}
                    continue

                # (c) L > 1  →  reuse or build a chain of length (L-1)
                if L in chain_cache:
                    # we have already built this chain once, just re-use its tail
                    tail = chain_cache[L][-1] if chain_cache[L] else u
                    G_unit.setdefault(tail, {})[v] = {1: {1}}
                    continue

                # otherwise build it now and store in cache
                prev = u
                helpers = []
                for _ in range(L - 1):
                    H = next_latent; next_latent += 1
                    G_unit.setdefault(prev, {})[H] = {1: {1}}
                    helpers.append(H)
                    prev = H
                # helpers[-1] is the tail of the chain
                chain_cache[L] = helpers
                G_unit.setdefault(prev, {})[v] = {1: {1}}

    # ── Step 3: ensure every observed node still appears ───────────────
    for u in obs_graph:
        G_unit.setdefault(u, {})

    return G_unit


def find_bcliques(graph):
    """
    Identify all B-cliques in the graph.
    A B-clique is any pair of non-empty subsets (V1, V2) such that
    every v1 in V1 has at least one non-empty edge (directed or bidirected)
    to every v2 in V2.
    """
    nodes = sorted(
        set(graph.keys()) |
        {nbr for nbrs in graph.values() for nbr in nbrs.keys()}
    )
    bcliques = []
    # iterate over all non-empty subsets for parents (V1)
    for r in range(1, len(nodes) + 1):
        for subset1 in combinations(nodes, r):
            V1 = set(subset1)
            # candidates for children: those v2 each v1 in V1 points to with a non-empty lag set
            candidates = [
                v2 for v2 in nodes
                if all(
                    v2 in graph.get(v1, {}) and any(graph[v1][v2].values())
                    for v1 in V1
                )
            ]
            # now choose any non-empty subset of those candidates as V2
            for s in range(1, len(candidates) + 1):
                for subset2 in combinations(candidates, s):
                    V2 = set(subset2)
                    bcliques.append((V1, V2))
    return bcliques

def find_maximal_bcliques(graph):
    """
    Identify maximal B-cliques in the graph.
    :param graph: Input graph (dictionary representation).
    :return: List of maximal B-cliques [(V1, V2)].
    """
    all_bcliques = find_bcliques(graph)  # Find all B-cliques
    maximal_bcliques = []

    for bc1 in all_bcliques:
        V1_1, V2_1 = bc1
        is_maximal = True

        for bc2 in all_bcliques:
            if bc1 == bc2:
                continue
            V1_2, V2_2 = bc2
            # Check if bc1 is a subset of bc2
            if V1_1.issubset(V1_2) and V2_1.issubset(V2_2):
                is_maximal = False
                break

        if is_maximal:
            maximal_bcliques.append(bc1)

    return maximal_bcliques

def compute_edge_lags(graph, bclique):
    V1, V2 = bclique
    edge_lags = {}
    for v1 in V1:
        for v2 in V2:
            ed = graph[v1][v2]
            # collect all lags under every etype present
            lags = set()
            for etype, raw in ed.items():
                if isinstance(raw, set):
                    lags |= raw
                else:
                    forest = raw if isinstance(raw, list) else [raw]
                    lags |= forest_to_set(forest)
            edge_lags[(v1, v2)] = lags
    return edge_lags

def compute_directed_lags(g, v, w, unobserved):
    """
    Compute d_O[v, w] using CRD paths.
    """
    paths = find_crd_paths(g, v, w, unobserved)  # Implement pathfinding logic
    edge_lags = set()
    for path in paths:
        edge_lags.update(compute_path_length(path))  # Aggregate path lengths
    return edge_lags

def compute_bi_directed_lags(g, v, w, unobserved):
    """
    Compute b_O[v, w] using h-treks.
    """
    treks = find_h_treks(g, v, w, unobserved)  # Implement h-trek logic
    edge_lags = set()
    for trek in treks:
        edge_lags.update(compute_trek_lags(trek))  # Aggregate trek lags
    return edge_lags

class SolutionNotFoundInTime(Exception):
    pass


def ptloopnum(pt):
    """
    Given a PathTree object returns the number of loops in it
    :param pt: PathTree object
    :return: number of loops (n)
    """

    def ptn(pt, n=0):
        for e in pt.loopset:
            if type(e) is int:
                n += 1
                continue
            n += ptn(e, n=1)
        return n

    return ptn(pt)


def ptnodenum(pt):
    """
    Given a PathTree object returns the number of latents that comprise it
    :param pt: PathTree object
    :return: number of nodes (n)
    """
    n = pt.preset - 1

    def ptn(pt, n=0):
        for e in pt.loopset:
            if type(e) is int:
                n += e - 1
                continue
            n += ptn(e, n=1)
        return n

    return n + ptn(pt)


def ptelement(pt, w):
    """
    An element generated by a PathTree with a given weight setting.
    :param pt: PathTree
    :param w: a list of weight IntVars (or nested lists thereof)
    :return: an IntExpr (or int, if there are no loops)
    """
    # 1) unwrap the root preset into a plain int
    if isinstance(pt.preset, set):
        if len(pt.preset) != 1:
            raise ValueError(f"ptelement: ambiguous preset {pt.preset}")
        base = next(iter(pt.preset))
    else:
        base = pt.preset

    def sumloops(node, weights):
        acc = 0
        children = list(node.loopset)
        for i, child in enumerate(children):
            if isinstance(child, int):
                # a flat loop
                acc += weights[i] * child
            else:
                # nested PathTree: unwrap its preset too
                if isinstance(child.preset, set):
                    assert len(child.preset) == 1
                    child_base = next(iter(child.preset))
                else:
                    child_base = child.preset

                # weights[i] is itself a pair [w_i, subweights]
                w_i, sub_w = weights[i]
                # multiply top‐level activation by child_base
                acc += w_i * child_base
                # plus “min(1,w_i) * recursive sum”
                acc += min(1, w_i) * sumloops(child, sub_w)
        return acc

    # if there are no loops, we just return base
    if not pt.loopset:
        return base
    # otherwise we add in all the loop contributions
    return base + sumloops(pt, w)


def weights_pt(pt, weights):
    c = [0]

    def crawl(pt, w, c):
        wl = []
        for e in pt.loopset:
            if type(e) is int:
                wl.append(w[c[0]])
                c[0] += 1
                continue
            ww = w[c[0]]
            c[0] += 1
            wl.append([ww, crawl(e, w, c)])
        return wl

    return crawl(pt, weights, c)


def extraloops_pt(pt, loops):  # loops are tuples (loop, weight)
    c = [0]

    def crawl(pt, l, c):
        first = [l[c[0]]]
        wl = []
        for e in pt.loopset:
            c[0] += 1
            if type(e) is int:
                wl.append(l[c[0]])
                continue
            wl.append(crawl(e, l, c))
        return first + [wl]

    return crawl(pt, loops, c)

def ptelement_extraloop(pt, w, eloops, solver):
    """
    Like ptelement, but adds in one extra “top-level” self-loop (eloops[0])
    and then the usual children/sub-loops (eloops[1]). Returns an IntExpr.
    """
    # 1) unwrap root preset into an int
    if isinstance(pt.preset, set):
        if len(pt.preset) != 1:
            raise ValueError(f"ambiguous preset {pt.preset}")
        base = next(iter(pt.preset))
    else:
        base = pt.preset

    # 2) Unpack the top-level extra-loop tuple: (L0, w0)
    L0, w0 = eloops[0]
    expr = base + w0 * L0

    def sumloops(node, weights, loopspecs):
        acc = 0
        children = list(node.loopset)
        child_specs = loopspecs[1]  # list of (Li, wi_flag) or [(Li, wi_flag), sub_specs]

        for i, child in enumerate(children):
            spec = child_specs[i]
            if isinstance(spec, tuple):
                Li, wi_flag = spec
                sub_specs = None
            else:
                (Li, wi_flag), sub_specs = spec

            if isinstance(child, int):
                # flat loop: w_i * child + w_i_flag * Li
                acc += weights[i] * child
                acc += wi_flag * Li
            else:
                # nested PathTree
                if isinstance(child.preset, set):
                    assert len(child.preset) == 1
                    cbase = next(iter(child.preset))
                else:
                    cbase = child.preset

                w_top, sub_w = weights[i]
                # top-level multiply + nested extra-loop
                acc += w_top * cbase
                nested = sumloops(child, sub_w, sub_specs) if sub_specs else 0
                acc += wi_flag * (Li + nested)

        return acc

    if pt.loopset:
        expr = expr + sumloops(pt, w, eloops)

    return expr



def isptelement_el(el, pt, w, eloops):
    return el == ptelement_extraloop(pt, w, eloops)


def isptsubset_el(elist, pt, w, eloops):
    for i in range(elist[-1]):
        if isptelement_el(i, pt, w, eloops):
            if not i in elist:
                return False
    return True


def isrightpt(el, elist, pt, w, eloops, solver):
    for i in range(elist[-1]):
        if isptelement_el(i, pt, w, eloops):
            if not i in elist:
                return False
        if i == el and not isptelement_el(i, pt, w, eloops):
            return False
    return True

def safe_int(x):
    """
    Attempt to convert x to an integer. If x is not a plain number (e.g. an InfiniteExpression),
    then raise an error.
    """
    if isinstance(x, (int, float)):
        return int(x)
    else:
        # We no longer attempt to convert symbolic types via int.
        raise TypeError(f"Cannot convert {x} (type {type(x)}) to int.")

def ptelements(pt, observed, seqlen=100, verbose=False, maxloop=100):
    """
    Generate the first `seqlen` delay values from a PathTree,
    but only those delays that are in the observed set.
    
    :param pt: A PathTree object.
    :param observed: A set (or list) of allowed delay values.
    :param seqlen: The number of delay values to generate.
    :param verbose: If True, prints solver statistics.
    :param maxloop: Upper bound for weight variables.
    :return: A list of delay values (only those that are in observed).
    """
    solver = pywrapcp.Solver("pt-elements")
    
    # If there are no children, the overall delay is just the preset.
    if not pt.loopset:
        # Only treat as "flat" if the preset is already a plain number.
        if isinstance(pt.preset, (int, float)):
            if pt.preset in observed:
                return [pt.preset]
            else:
                raise ValueError(f"Preset {pt.preset} is not in the observed set {observed}.")
        # Otherwise, fall through to constraint solving.
    
    # Otherwise, use the constraint solver.
    N = ptloopnum(pt)  # Number of loops in the PathTree
    weights = [solver.IntVar(0, maxloop, f"w[{i:04d}]") for i in range(N)]
    #weights = [solver.IntVar(0, maxloop) for i in range(N)]
    ws = weights_pt(pt, weights)
    delay_expr = ptelement(pt, ws)
    
    # Make sure the observed set contains plain numbers.
    numeric_observed = [int(x) for x in observed if isinstance(x, (int, float))]
    if not numeric_observed:
        raise ValueError("Observed set does not contain any numeric values.")
    
    # Add a constraint that the computed delay must be one of the observed values.
    solver.Add(solver.MemberCt(delay_expr, numeric_observed))
    
    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)
    
    num_solutions = 0
    els = set()
    while solver.NextSolution():
        w_values = [x.Value() for x in weights]
        # weights_pt can work with a plain list of numbers
        nested_ws = weights_pt(pt, w_values)
        delay_val = ptelement(pt, nested_ws)
        els.add(delay_val)
        num_solutions += 1
        if len(els) >= seqlen:
            break
    solver.EndSearch()
    
    if verbose:
        print("num_solutions:", num_solutions)
        print("failures:", solver.Failures())
        print("branches:", solver.Branches())
        print("WallTime:", solver.WallTime())
    
    return list(els)


def isptelement(pt, element, verbose=False, maxloop=100):
    """
    Check if an integer element is in the weight set represented by the PathTree.
    If the PathTree is flat (no children), and its preset is a set, check membership;
    otherwise, use the constraint solver to verify if 'element' is generated by the PathTree.
    
    :param pt: a PathTree object
    :param element: an integer candidate delay
    :param verbose: whether to print debugging information
    :param maxloop: upper bound for weight variables
    :return: True if 'element' is representable by pt, else False
    """
    solver = pywrapcp.Solver("isptelement")
    weights = []
    N = ptloopnum(pt)  # Number of loops (children) in the PathTree

    # If there are no children, then the delay is given directly by pt.preset.
    if not N:
        if isinstance(pt.preset, set):
            return element in pt.preset
        else:
            return element == pt.preset

    # Otherwise, initialize weight variables.
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))
        #weights.append(solver.IntVar(0, maxloop))

    # Get the weighted representation of the PathTree.
    wpt = weights_pt(pt, weights)
    solver.Add(element == ptelement(pt, wpt))

    solution = solver.Assignment()
    solution.Add(weights)
    db = solver.Phase(weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    solution_exists = False
    while solver.NextSolution():
        solution_exists = True
        break
    solver.EndSearch()

    if verbose:
        print("Failures:", solver.Failures())
        print("Branches:", solver.Branches())
        print("WallTime:", solver.WallTime())

    return solution_exists


def loops_and_weights(solver, loops, weights):
    """
    Add constraints to solver that make sure loops are not generated if subtree is not active due to a zero weight upstream
    :param solver:
    :param loops:
    :param weights:
    :return:
    """

    def recurse(s, l, w):
        for ww, ll in zip(w, l):
            if type(ww) is list:
                for e in flatten(ll):
                    s.Add((ww[0] == 0) <= (e == 0))
                recurse(s, ll[1:], ww[1:])
            else:
                for e in flatten(ll):
                    s.Add((ww == 0) <= (e == 0))

    recurse(solver, loops[1], weights)


def eloops_simplify(eloops):
    l = []
    for e in eloops:
        if type(e) is list:
            l.append(eloops_simplify(e))
        else:
            l.append(int(e[0].Value()))
    return l

def ptaugmented(pt, eloops):
    """
    Given a PathTree pt and solver-defined loopspecs eloops:
      eloops[0] = (L0, w0)        # top-level loop var
      eloops[1] = [(L1, w1), ...] # one pair per original loop in pt.loopset
    Return a new PathTree where:
      - if w0>0, the new loop L0 is added flat on the root
      - for each original loop Li, if wi>0, that Li is *nested* under its original loop
    """

    def unwrap(x):
        # pull plain int out of IntVar or singleton set
        if hasattr(x, "Value"):
            return x.Value()
        if isinstance(x, set):
            return next(iter(x))
        return x

    base = unwrap(pt.preset)

    def augment(node, loopspecs):
        new_loops = set()

        # 1) possibly add a brand-new loop on the root
        L0, w0 = loopspecs[0]
        if unwrap(w0):
            new_loops.add(unwrap(L0))

        # 2) for each original child‐loop, either keep it flat or nest the extra loop under it
        for (orig_loop, (Li, wi)) in zip(node.loopset, loopspecs[1]):
            orig_len = unwrap(orig_loop if isinstance(orig_loop, int) else orig_loop.preset)
            if unwrap(wi):
                # nest the Li under that orig_loop
                new_loops.add(
                    PathTree(preset={orig_len}, loopset={unwrap(Li)})
                )
            else:
                # keep the old loop intact
                new_loops.add(orig_len)

        return PathTree(preset={base}, loopset=new_loops)

    return augment(pt, eloops)


def ptsubset(pt, elist):
    """
    Return True if every delay that the PathTree `pt` can generate—up to the
    maximum observed lag—is contained in the observed set `elist`.
    """
    # Determine the maximum lag to check
    max_lag = max(elist)
    # For each possible generated delay up to max_lag
    for i in range(max_lag + 1):
        # If pt can generate i but i is not in the observed set, reject
        if isptelement(pt, i) and i not in elist:
            return False
    return True


def smallest_pt(ptlist):
    if ptlist:
        idx = np.argsort(map(ptnodenum, ptlist))
        sol = ptlist[idx[0]]
    else:
        sol = None
    return sol


def pairprint(pt1, pt2, k=40):
    print(np.c_[pt2seq(pt1, k), pt2seq(pt2, k)])


def etesteq(pt1, pt2, k=100):
    a1 = np.asarray(pt2seq(pt1, k))
    a2 = np.asarray(pt2seq(pt2, k))
    return np.sum(a1 - a2) == 0


def keeptreegrow(pt, e, seq, cutoff=10, cap=1000):
    t = None
    while t is None:
        t = growtree(pt, e, seq, cutoff=cutoff)
        cutoff += 10
        if cutoff > cap:
            raise SolutionNotFoundInTime("Cannot keep the tree growing")
    return t


def add_element(d, pt):
    """
    Add a PathTree to dictionary d such that it is either appended to the list or added anew
    Args:
        d: a dictionary
        pt: a PathTree

    Returns:

    """
    key = ptnodenum(pt)
    if key in d:
        d[key].append(pt)
    else:
        d[key] = pt


def del_element(d, pt, key=None):
    """
    Delete a PathTree from dictionary d such that it is either removed from the list or the list that only contains one element is removed
    Args:
        d: a dictionary
        pt: a PathTree

    Returns:

    """
    if key is None:
        key = ptnodenum(pt)
    if len(d[key]) == 1:
        del d[key]
    else:
        d[key].remove(pt)


def swap_elements(d, pt1, pt2, key=None):
    del_element(d, pt1, key=key)
    add_element(d, pt2)


def seq2pt(seq, verbose=False, cutoff=100):
    if not seq:
        return None

    pt = PathTree({}, pre=seq[0])
    pts = SortedDict()  # PathTrees
    pts[ptnodenum(pt)] = [pt]

    for e in seq[1:]:
        e_is_in = False
        for key in pts:
            for pt in pts[key]:
                if verbose:
                    print(e)
                try:
                    newpt = keeptreegrow(pt, e, seq, cutoff=cutoff)
                    swap_elements(pts, pt, newpt, key=key)
                    e_is_in = True
                    break
                except SolutionNotFoundInTime:
                    continue
        if not e_is_in:
            newpt = PathTree({}, pre=e)
            add_element(d, newpt)

    return pt

def update_edge_lags(graph, bclique, new_pt):
    V1, V2 = bclique

    # 5) overwrite each edge in the B-clique to use this tree
    for u in V1:
        for v in V2:
            graph[u][v][1] = new_pt

    return graph

def growtree(pt, element, ref_elements, verbose=False, maxloop=100, cutoff=100):
    """
    Add a loop with the minimal length to a path tree to enable it to generate a given element and still be a subset of a given list
    :param pt: a path tree object from pathtree.py
    :param element: an integer to check for presence in the weight
    :param ref_elements: a (finite) list that should be a superset of numbers generated by the new path tree, for numbers smaller than tosubset[-1]
    :param verbose: whether to print debugging information
    :return: a PathTree augmented with a new loop
    """
    solver = pywrapcp.Solver("loop_an_element")

    # PathTree already can generate that number. Just to foolproof
    if isptelement(pt, element):
        return pt

    # declare variables
    weights = []  # weights denoting how many times a loop is active (marginalized)
    loops = []  # extra loops that can be potentially added
    lweights = []  # weights for the extra loops (marginalized out in the end)
    ltuples = []  # tuple list to hold loops and weights together

    N = ptloopnum(pt)  # number of loops in the PathTree
    for i in range(N):
        weights.append(solver.IntVar(0, maxloop, "w[%04i]" % i))

    for i in range(N + 1):
        w = solver.IntVar(0, maxloop, "lw[%04i]" % i)
        l = solver.IntVar(0, maxloop, "l[%04i]" % i)
        lweights.append(w)  # loop related weight
        loops.append(l)
        ltuples.append((l, w))

    eloops = extraloops_pt(pt, ltuples)
    ws = weights_pt(pt, weights)

    # declare constraints
    # 1) produce the IntExpr
    expr = ptelement_extraloop(pt, ws, eloops, solver)
    # 2) cast your allowed‐lags set into a plain list of ints
    allowed = list(ref_elements)           # e.g. [6,8,10,11,12,13]
    # 3) add the membership constraint
    solver.Add(solver.MemberCt(expr, allowed))
    solver.Add(element == ptelement_extraloop(pt, ws, eloops, solver))  # make sure the element can be generated
    solver.Add(solver.Count(loops, 0, len(loops) - 1))  # only one loop is on
    solver.Add(solver.Count(lweights, 0, len(lweights) - 1))  # only one loop is weighted
    for i in range(len(lweights)):
        solver.Add((lweights[i] == 0) <= (loops[i] == 0))  # if a loop has weight zero then it can't be active
        # solver.Add(lweights[i] >= loops[i])
    loops_and_weights(solver, eloops, ws)  # if a subtree is off (weight zero) no need to add loops

    # run the solver
    solution = solver.Assignment()
    solution.Add(loops)
    db = solver.Phase(loops + lweights + weights,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)
    solver.NewSearch(db)

    numsol = 0
    pts = []
    while solver.NextSolution():
        # print numsol,
        new_pt = ptaugmented(pt, eloops)
        if verbose:
            print("trying PathTree: ", new_pt)
        if ptsubset(new_pt, ref_elements):
            pts.append(new_pt)
            if verbose:
                print("OK PathTree: ", pts[-1])
        numsol += 1
        if numsol >= cutoff:
            break
    solver.EndSearch()

    # output solutions
    if verbose:
        print("solutions:", numsol)
        print("failures:", solver.Failures())
        print("branches:", solver.Branches())
        print("WallTime:", solver.WallTime())
        print("for ", element, "solutions found ", numsol)

    return smallest_pt(pts)

def apply_minimal_refinement(graph, bcliques, cap=5):
    for V1, V2 in bcliques:
        for v1 in V1:
            for v2 in V2:
                if v2 not in graph[v1]:
                    continue
                for etype, raw in list(graph[v1][v2].items()):

                    # ── NEW: unwrap the common ‘flat preset = {…}’ pattern ──
                    if (isinstance(raw, PathTree) and not raw.loopset
                        and isinstance(raw.preset, set) and len(raw.preset) > 1):
                        lag_set = raw.preset
                    else:
                        if isinstance(raw, set):
                            lag_set = raw
                        else:                      # PathTree or list thereof
                            forest = raw if isinstance(raw, list) else [raw]
                            lag_set = forest_to_set(forest, cap=cap)
                    # ----------------------------------------------------------------

                    if not lag_set:
                        continue

                    a = min(lag_set)
                    diffs = [x - a for x in lag_set]
                    step  = gcd(*diffs) if len(diffs) > 1 else 0

                    # arithmetic progression → single-loop *only* if it stays inside lag_set
                    if step and lag_set == {a + k*step for k in range(len(lag_set))}:
                        cand = PathTree(preset={a}); cand.add_loop(step)
                        gen  = forest_to_set(cand, cap=cap)
                        if gen.issubset(lag_set):          # ✓ safe -- keep it
                            graph[v1][v2][etype] = cand
                            continue                       # done
                    # else fall through to “explicit” representation

                    # original heuristics
                    if len(lag_set) == 2:
                        b = max(lag_set)
                        # candidate PathTree  (a,<b-a>)
                        pt = PathTree(preset={a}); pt.add_loop(b - a)

                        # keep it *only* if it generates no new lags
                        gen = forest_to_set(pt, cap=cap)      # e.g. {a, a+d, a+2d, …}
                        if gen.issubset(lag_set):
                            graph[v1][v2][etype] = pt
                        else:          # otherwise fall back to “one tree per lag”
                            graph[v1][v2][etype] = [
                                PathTree(preset={e}) for e in sorted(lag_set)
                            ]
                    else:
                        graph[v1][v2][etype] = [
                            PathTree(preset={e}) for e in sorted(lag_set)
                        ]
    return graph

def pt2seq(pt, num):
    if not pt.children:
        return [pt.preset]
    i = 0
    s = set()
    while len(s) < num:
        if isptelement(pt, i, maxloop=10 * num):
            s.add(i)
        i += 1
    l = list(s)
    l.sort()
    return l


def s2spt(s):  # convert edge set to pt
    ss = set()
    for e in s:
        if type(e) is int:
            ss.add(PathTree({0}, pre={e}))
            continue
        ss.add(e)
    return ss


def spt_elements(spt, num):
    """
    Generate numbers from a set of PathTrees
    :param spt: set of PathTrees
    :param num: number of elements (from the first) to generate
    :return: list of num numbers
    """
    i = 0
    s = set()
    while len(s) < num:
        if issptelement(spt, i):
            s.add(i)
        i += 1
    return list(s)


def issptelement(spt, element):
    a = False
    for pt in s2spt(spt):
        a = a or isptelement(pt, element)
    return a

def unify_edge_representation(graph):
    """
    Convert all edges in the graph so that they store a PathTree, even if
    originally they were just a set of integers.
    """
    for u in graph:
        for v in graph[u]:
            edge_data = graph[u][v]
            if 1 in edge_data:  # directed
                if isinstance(edge_data[1], set):
                    # Wrap the set in a PathTree
                    edge_data[1] = PathTree(preset=edge_data[1], loopset=set())
            elif 2 in edge_data:  # bidirected
                if isinstance(edge_data[2], set):
                    edge_data[2] = PathTree(preset=edge_data[2], loopset=set())
    return graph

def print_graph_as_pathtrees(graph):
    """
    Print each edge in the graph on a separate line.
    For each edge from u to v, if a directed edge (key 1) exists,
    print "u -> v = ..." on its own line.
    If a bidirected edge (key 2) exists, print "u <-> v = ..." on a separate line.
    """
    for u in sorted(graph.keys()):
        for v in sorted(graph[u].keys()):
            edge_data = graph[u][v]
            if 1 in edge_data:
                print(f"{u} -> {v} = {edge_data[1]}")
            if 2 in edge_data:
                print(f"{u} <-> {v} = {edge_data[2]}")
