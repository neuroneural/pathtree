#!/usr/bin/env bash
set -u
set -o pipefail

INPUT_DIR="ringmore_inputs"
GSET_DIR="ringmore_gsets"
OUT_DIR="ringmore_reverse"
SUMMARY_FILE="$OUT_DIR/summary.txt"

CORE="hide_nodes_program.lp"
PRUNE_REV="prune_reverse.lp"
OUT_REV="output_directives.lp"     # ok even if it shows out_dir/out_bi too

# Optional forward check (set to 1 to enable)
DO_FORWARD_CHECK=1
PRUNE_FWD="prune_forward.lp"
OUT_FWD="output_directives.lp"

mkdir -p "$OUT_DIR"
true_count=0
false_count=0
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
  # canonicalize a file of atoms:
  # - one atom per line
  # - strip blanks
  # - sort + unique
  tr ' ' '\n' | sed '/^$/d' | sort -u
}

echo "Running reverse on ringmore lagsets"
echo "Input graphs: $INPUT_DIR"
echo "Gsets:        $GSET_DIR"
echo "Output dir:   $OUT_DIR"
echo

for gset in "$GSET_DIR"/*_gset.lp; do
  [ -f "$gset" ] || continue

  base=$(basename "$gset" _gset.lp)              # input_001
  orig="$INPUT_DIR/${base}.lp"                   # ringmore_inputs/input_001.lp

  inv="$OUT_DIR/${base}_inverse_input.lp"
  raw="$OUT_DIR/${base}_candidate_raw.lp"
  raw_lines="$OUT_DIR/${base}_candidate_raw_lines.lp"
  cand_graph="$OUT_DIR/${base}_candidate_graph.lp"

  err="$OUT_DIR/${base}.err"

  gset2="$OUT_DIR/${base}_gset2.lp"
  same="$OUT_DIR/${base}_same.txt"
  fwd_err="$OUT_DIR/${base}_forward2.err"

  echo "→ $base"

  if [ ! -f "$orig" ]; then
    echo "  ERROR: missing original input file $orig"
    exit 1
  fi

  # --- build inverse input: observed + preserved observed->observed edges + lagset facts ---
  : > "$inv"

  grep '^observed(' "$orig" >> "$inv"

  # preserve ONLY observed->observed unit edges from original input
  awk '
    /^observed\(/ { obs[$0]=1 }
    /^edge\(/ {
      split($0,a,"[(),]")
      u=a[2]; v=a[3]
      if (obs["observed(" u ")."] && obs["observed(" v ")."])
        print $0
    }
  ' "$orig" >> "$inv"

  # target lagset facts go in verbatim (out_dir/out_bi)
  cat "$gset" >> "$inv"

  # --- reverse solve ---
  set +e
  clingo "$inv" "$CORE" "$PRUNE_REV" "$OUT_REV" --outf=1 --opt-mode=optN --models=1 2>"$err" \
    | atoms_only > "$raw"
  rc=$?
  set -e

  # rc: 10 SAT, 20 UNSAT, 30 UNKNOWN
  if [ "$rc" -ne 10 ] && [ "$rc" -ne 20 ] && [ "$rc" -ne 30 ]; then
    echo "  ERROR: clingo failed for $base (rc=$rc). See $err"
    exit 1
  fi
  if [ "$rc" -eq 20 ]; then
    echo "  WARNING: UNSAT for $base"
  elif [ "$rc" -eq 30 ]; then
    echo "  WARNING: UNKNOWN for $base"
  fi

  # normalize raw to one atom per line (handles single-line outputs)
  tr ' ' '\n' < "$raw" | sed '/^$/d' > "$raw_lines"

  # --- build candidate_graph in the old style: observed/hidden/edge only ---
  : > "$cand_graph"

  # observed facts from original
  grep '^observed(' "$orig" >> "$cand_graph"

  # preserve ONLY observed->observed edges from original
  awk '
    /^observed\(/ { obs[$0]=1 }
    /^edge\(/ {
      split($0,a,"[(),]")
      u=a[2]; v=a[3]
      if (obs["observed(" u ")."] && obs["observed(" v ")."])
        print $0
    }
  ' "$orig" >> "$cand_graph"

  # hidden facts from out_active(...)
  awk '
    /^out_active\(/ {
      line=$0
      sub(/^out_active\(/,"hidden(",line)
      print line
    }
  ' "$raw_lines" >> "$cand_graph"

  # edges from out_edge(...)
  awk '
    /^out_edge\(/ {
      line=$0
      sub(/^out_edge\(/,"edge(",line)
      print line
    }
  ' "$raw_lines" >> "$cand_graph"

  # --- optional forward check: hide candidate_graph and compare to original gset ---
  if [ "$DO_FORWARD_CHECK" -eq 1 ]; then
    set +e
    clingo "$cand_graph" "$CORE" "$PRUNE_FWD" "$OUT_FWD" --outf=1 2>"$fwd_err" \
      | atoms_only > "$gset2"
    fwd_rc=$?
    set -e

    if [ "$fwd_rc" -ne 10 ] && [ "$fwd_rc" -ne 20 ] && [ "$fwd_rc" -ne 30 ]; then
      echo "  ERROR: forward-check clingo failed for $base (rc=$fwd_rc). See $fwd_err"
      exit 1
    fi

    canon1="$OUT_DIR/${base}_gset1.canon"
    canon2="$OUT_DIR/${base}_gset2.canon"

    canon_lp < "$gset"  > "$canon1"
    canon_lp < "$gset2" > "$canon2"

    if diff -q "$canon1" "$canon2" >/dev/null; then
      echo True > "$same"
      true_count=$((true_count+1))
    else
      echo False > "$same"
      false_count=$((false_count+1))
    fi

  fi

  echo "  wrote: $(basename "$raw"), $(basename "$cand_graph")"
done

{
  echo "Done. Outputs are in $OUT_DIR"
  echo "True  : $true_count"
  echo "False : $false_count"
} | tee "$SUMMARY_FILE"
