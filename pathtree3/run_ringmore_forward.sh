#!/usr/bin/env bash
set -u
set -o pipefail

INPUT_DIR="ringmore_inputs"
OUT_DIR="ringmore_gsets"

CORE="hide_nodes_program.lp"
PRUNE_FWD="prune_forward.lp"
OUT_FWD="output_directives.lp"   # better: output_forward.lp with only out_dir/out_bi

mkdir -p "$OUT_DIR"

echo "Running forward hide on ringmore graphs"
echo "Input dir:  $INPUT_DIR"
echo "Output dir: $OUT_DIR"
echo

for f in "$INPUT_DIR"/input_*.lp; do
  base=$(basename "$f" .lp)
  out="$OUT_DIR/${base}_gset.lp"
  err="$OUT_DIR/${base}.err"

  echo "→ $base"

  set +e
  clingo "$f" "$CORE" "$PRUNE_FWD" "$OUT_FWD" --outf=1 2>"$err" \
    | awk '
        $0 ~ /^%/ { next }
        $0 == "ANSWER" { next }
        $0 ~ /^COST/ { next }
        $0 == "OPTIMUM" { next }
        NF == 0 { next }
        { print }
      ' > "$out"
  rc=$?
  set -e

  # clingo rc: 10 SAT, 20 UNSAT, 30 UNKNOWN
  if [ "$rc" -ne 10 ] && [ "$rc" -ne 20 ] && [ "$rc" -ne 30 ]; then
    echo "  ERROR: clingo failed for $base (rc=$rc). See $err"
    exit 1
  fi

  if [ "$rc" -eq 30 ]; then
    echo "  WARNING: clingo returned UNKNOWN (rc=30) for $base"
  fi

  nlines=$(wc -l < "$out" | tr -d " ")
  echo "  wrote $nlines lines (rc=$rc)"
done

echo
echo "Done. Lagsets written to $OUT_DIR"
