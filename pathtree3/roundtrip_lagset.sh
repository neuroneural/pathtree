#!/usr/bin/env bash
set -u

INPUT_LP="${1:-input.lp}"

HIDE_CORE="hide_nodes_program.lp"
PRUNE_FWD="prune_forward.lp"
PRUNE_REV="prune_reverse.lp"
OUT_LP="output_directives.lp"

atoms_only() {
  awk '
    $0 ~ /^%/ { next }
    $0 == "ANSWER" { next }
    $0 ~ /^COST/ { next }
    $0 == "OPTIMUM" { next }
    NF == 0 { next }
    { print }
  '
}
canon_lp() {
  tr ' ' '\n' | sed '/^$/d' | sort -u
}
die() {
  echo "ERROR: $1" >&2
  exit 1
}

echo "[0/6] Checking files"
[ -f "$INPUT_LP" ]   || die "Missing $INPUT_LP"
[ -f "$HIDE_CORE" ]  || die "Missing $HIDE_CORE"
[ -f "$PRUNE_FWD" ]  || die "Missing $PRUNE_FWD"
[ -f "$PRUNE_REV" ]  || die "Missing $PRUNE_REV"
[ -f "$OUT_LP" ]     || die "Missing $OUT_LP"
command -v clingo >/dev/null 2>&1 || die "clingo not found"

echo "[1/6] Forward hide: ${INPUT_LP} -> gset1.lp"
clingo "$INPUT_LP" "$HIDE_CORE" "$PRUNE_FWD" "$OUT_LP" --outf=1 2>clingo_forward1.err \
  | atoms_only > gset1.lp

if [ ! -s gset1.lp ]; then
  echo "Forward produced empty gset1.lp. Stderr was:" >&2
  sed -n '1,120p' clingo_forward1.err >&2
  die "gset1.lp is empty"
fi

echo "[2/6] Build reverse_facts.lp (observed + preserved obs->obs edges + target lagset facts)"
# Observed nodes
grep '^observed(' "$INPUT_LP" > reverse_facts.lp
[ -s reverse_facts.lp ] || die "No observed(...) facts found in $INPUT_LP"

# Preserve ONLY observed->observed unit edges from the original input
awk '
  /^observed\(/ { obs[$0]=1 }
  /^edge\(/ {
    split($0,a,"[(),]")
    u=a[2]; v=a[3]
    if (obs["observed(" u ")"] && obs["observed(" v ")"])
      print $0
  }
' "$INPUT_LP" >> reverse_facts.lp

# Target lagset = gset1 verbatim (already out_dir/out_bi)
cat gset1.lp >> reverse_facts.lp

echo "[3/6] Reverse solve -> candidate_raw.lp (one optimal model)"
clingo reverse_facts.lp "$HIDE_CORE" "$PRUNE_REV" "$OUT_LP" --opt-mode=optN --models=1 --outf=1 2>clingo_reverse.err \
  | atoms_only > candidate_raw.lp

if [ ! -s candidate_raw.lp ]; then
  echo "Reverse produced empty candidate_raw.lp. Stderr was:" >&2
  sed -n '1,200p' clingo_reverse.err >&2
  die "candidate_raw.lp is empty"
fi

echo "[4/6] Normalize candidate_raw.lp to one atom per line"
tr ' ' '\n' < candidate_raw.lp | sed '/^$/d' > candidate_raw_lines.lp

echo "[5/6] Build candidate_graph.lp (facts for second forward run)"
# Start with observed nodes and preserved observed->observed edges
grep -E '^(observed|edge)\(' reverse_facts.lp > candidate_graph.lp

# Add hidden nodes selected by reverse + all chosen edges
awk '
  /^out_active\(/ {
    line=$0; sub(/^out_active\(/,"hidden(",line); print line; next
  }
  /^out_edge\(/ {
    line=$0; sub(/^out_edge\(/,"edge(",line); print line; next
  }
' candidate_raw_lines.lp >> candidate_graph.lp

grep -q '^edge(' candidate_graph.lp || die "candidate_graph.lp has no edge(...) facts"

echo "[6/6] Forward hide candidate_graph.lp -> gset2.lp and compare"
clingo candidate_graph.lp "$HIDE_CORE" "$PRUNE_FWD" "$OUT_LP" --outf=1 2>clingo_forward2.err \
  | atoms_only > gset2.lp

if [ ! -s gset2.lp ]; then
  echo "Second forward produced empty gset2.lp. Stderr was:" >&2
  sed -n '1,120p' clingo_forward2.err >&2
  die "gset2.lp is empty"
fi

canon_lp < gset1.lp > gset1.canon
canon_lp < gset2.lp > gset2.canon

if diff -q gset1.canon gset2.canon >/dev/null; then
  echo True > same.txt
  echo "same.txt: True"
else
  echo False > same.txt
  echo "same.txt: False"
fi


rm -f candidate_raw_lines.lp
echo "Wrote: gset1.lp reverse_facts.lp candidate_raw.lp candidate_graph.lp gset2.lp same.txt"
echo "Logs: clingo_forward1.err clingo_reverse.err clingo_forward2.err"
