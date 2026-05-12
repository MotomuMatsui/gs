#!/usr/bin/env bash
# Regression check for gs2 against committed baseline.
#
# Two modes:
#   strict   (default) — OMP_NUM_THREADS=1, exact byte diff against baseline.
#                        Use for FP-preserving refactors.
#   tolerant           — default OMP, topology must match + EP support within --eps.
#                        Use when changes legitimately reorder FP accumulation
#                        (alternative eigensolvers, OMP scheduling changes).
#
# Picking --eps for tolerant mode: at ep_num=100 the binomial noise on each
# clade's EP is roughly 1/sqrt(100)=0.10, with occasional 2σ outliers up to
# ~0.20.  Use 0.20 for routine checks, 0.10 only with ep_num>=400.
#
# Usage:
#   bench/regression.sh                          # strict, bit-exact
#   bench/regression.sh --tolerant               # tolerant, eps=0.20
#   bench/regression.sh --tolerant --eps 0.10    # tighter tolerance
#   GS2=./gs2_gpu bench/regression.sh            # test the GPU build

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GS2="${GS2:-$ROOT/gs2}"
INPUT="$ROOT/example/200.faa"
BASELINE="$ROOT/bench/baseline/200_e100_seed12345.nwk"
COMPARE="$ROOT/bench/compare_newick.py"

MODE="strict"
EPS="0.20"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --strict)   MODE="strict"; shift ;;
    --tolerant) MODE="tolerant"; shift ;;
    --eps)      EPS="$2"; shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 64 ;;
  esac
done

if [[ ! -x "$GS2" ]]; then
  echo "ERROR: $GS2 not found or not executable. Build first: make" >&2
  exit 64
fi
if [[ ! -f "$BASELINE" ]]; then
  echo "ERROR: baseline missing: $BASELINE" >&2
  exit 64
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

OUT="$TMP/run.nwk"

echo "[regression] mode=$MODE  gs2=$GS2  input=$INPUT  baseline=$BASELINE"

if [[ "$MODE" == "strict" ]]; then
  # OPENBLAS_NUM_THREADS=1 keeps strict mode deterministic when gs2 is linked
  # against multi-threaded OpenBLAS (OpenBLAS's parallel reductions are not
  # FP-associative).  Harmless when linked against single-threaded BLAS.
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "$GS2" -s -e 100 -r 12345 "$INPUT" > "$OUT"
  if diff -q "$BASELINE" "$OUT" > /dev/null; then
    echo "OK: byte-identical to baseline (strict, OMP=1)"
    exit 0
  fi
  echo "FAIL: output differs from baseline."
  echo "  baseline: $BASELINE"
  echo "  actual:   $OUT  (kept for inspection: copying to /tmp/gs2_regression_fail.nwk)"
  cp "$OUT" /tmp/gs2_regression_fail.nwk
  diff "$BASELINE" "$OUT" | head -20
  exit 1
fi

# Tolerant mode: run with whatever OMP setting is in the environment.
"$GS2" -s -e 100 -r 12345 "$INPUT" > "$OUT"
python3 "$COMPARE" "$BASELINE" "$OUT" --eps "$EPS"
