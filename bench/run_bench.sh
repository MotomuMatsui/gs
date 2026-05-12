#!/usr/bin/env bash
# Wall-clock benchmark for gs2.
#
# Runs `./gs2 -s -e $E -r 12345 example/200.faa` for several EP iteration
# counts and reports the median wall-time over $REPEATS runs.
#
# Knobs (env vars):
#   GS2       path to the gs2 binary (default: ./gs2)
#   INPUT     FASTA input (default: example/200.faa)
#   E_VALUES  space-separated -e values (default: "0 10 50 100")
#   REPEATS   per-config repetition count (default: 3)
#   THREADS   value for OMP_NUM_THREADS (default: unset, gs2 picks)
#
# `-e 0` measures the GS-only path (no EP), useful for isolating spectral
# clustering cost from EP cost.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GS2="${GS2:-$ROOT/gs2}"
INPUT="${INPUT:-$ROOT/example/200.faa}"
E_VALUES="${E_VALUES:-0 10 50 100}"
REPEATS="${REPEATS:-3}"

if [[ ! -x "$GS2" ]]; then
  echo "ERROR: $GS2 not found or not executable. Build first: make" >&2
  exit 64
fi

if [[ -n "${THREADS:-}" ]]; then
  export OMP_NUM_THREADS="$THREADS"
fi

echo "# gs2 benchmark"
echo "# binary : $GS2"
echo "# input  : $INPUT"
echo "# omp    : ${OMP_NUM_THREADS:-(default)}"
echo "# repeats: $REPEATS"
echo "#"
echo "# NOTE: end-to-end wall time INCLUDES MMseqs2 / BLASTN alignment."
echo "#       For N=200 protein input this is ~95-97% of wall time."
echo "#       To isolate gs2 internals, materialize a similarity matrix"
echo "#       (see bench/README.md) and benchmark with 'gs2 -w simmat.txt'."
echo
printf "%-6s %-10s %-10s %-10s %s\n" "ep" "median_s" "min_s" "max_s" "runs"
printf "%-6s %-10s %-10s %-10s %s\n" "------" "----------" "----------" "----------" "----------"

# Use bash 5+ $EPOCHREALTIME (microsecond-resolution wallclock, no /usr/bin/time fork).
time_run() {
  local t0 t1
  t0="$EPOCHREALTIME"
  "$GS2" -s -e "$1" -r 12345 "$INPUT" > /dev/null
  t1="$EPOCHREALTIME"
  awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.4f\n", b - a}'
}

for E in $E_VALUES; do
  samples=()
  for ((r=1; r<=REPEATS; r++)); do
    samples+=("$(time_run "$E")")
  done

  # Sort numerically to get min / median / max.
  sorted=$(printf '%s\n' "${samples[@]}" | sort -n)
  mapfile -t arr <<< "$sorted"
  n=${#arr[@]}
  mid=$(( n / 2 ))
  if (( n % 2 == 0 )); then
    median=$(awk -v a="${arr[mid-1]}" -v b="${arr[mid]}" 'BEGIN{printf "%.4f", (a+b)/2}')
  else
    median="${arr[mid]}"
  fi
  min="${arr[0]}"
  max="${arr[n-1]}"

  printf "%-6s %-10s %-10s %-10s %s\n" "$E" "$median" "$min" "$max" "${samples[*]}"
done
