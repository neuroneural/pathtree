#!/usr/bin/env python3
import random
import argparse
from pathlib import Path

def parse_hidden(s: str):
    s = s.strip()
    if not s:
        return set()
    s = s.strip("{}")
    return {int(x.strip()) for x in s.split(",") if x.strip()}

def ring_edges(n: int):
    return {(i, i+1) for i in range(1, n)} | {(n, 1)}

def all_possible_edges(n: int):
    # self-loops allowed
    return {(u, v) for u in range(1, n+1) for v in range(1, n+1)}

def write_lp(path: Path, n: int, hidden: set[int], edges: set[tuple[int,int]]):
    lines = []
    for i in range(1, n+1):
        lines.append(f"{'hidden' if i in hidden else 'observed'}({i}).")
    lines.append("")
    for (u, v) in sorted(edges):
        lines.append(f"edge({u},{v},1).")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Ngraph", type=int, required=True)
    ap.add_argument("--Nnode", type=int, required=True)
    ap.add_argument("--hidden", type=str, default="")
    ap.add_argument("--Nextra", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", type=str, default="ringmore_inputs")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hidden = parse_hidden(args.hidden)
    n = args.Nnode
    base = ring_edges(n)

    candidate_pool = sorted(all_possible_edges(n) - base)

    for g in range(1, args.Ngraph + 1):
        edges = set(base)
        if args.Nextra > 0:
            if args.Nextra > len(candidate_pool):
                raise ValueError("Nextra is too large for the available non-ring edges.")
            extra = rng.sample(candidate_pool, args.Nextra)
            edges.update(extra)

        write_lp(outdir / f"input_{g:03d}.lp", n, hidden, edges)

if __name__ == "__main__":
    main()