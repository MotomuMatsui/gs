#!/usr/bin/env bash
# Profile gs2.
#
# Subcommands:
#   stat    — `perf stat -d` summary (cycles, IPC, cache, etc.) when permitted,
#             else falls back to `/usr/bin/time -v` (wall + RSS + page faults).
#   record  — `perf record -g`, then `perf report` interactive
#   flame   — perf record + folded stacks + FlameGraph SVG  (requires
#             FlameGraph scripts on PATH or FLAMEGRAPH_DIR set)
#
# Knobs (env vars):
#   GS2       path to gs2 (default: ./gs2)
#   INPUT     FASTA input (default: example/200.faa)
#   E         -e value (default: 100)
#   OUT       output file base name (default: bench/profile-<subcmd>)
#   FREQ      perf sampling frequency for record/flame (default: 999)
#   FLAMEGRAPH_DIR  path to FlameGraph scripts (stackcollapse-perf.pl, flamegraph.pl)
#
# perf permissions
# ----------------
# Most distros default to /proc/sys/kernel/perf_event_paranoid = 4, which
# blocks every perf-event for non-privileged users.  To use `stat`/`record`/
# `flame`, lower it (transient — reset on reboot):
#     sudo sysctl kernel.perf_event_paranoid=1
# Or, less invasive: run this script under sudo.
#
# Output files are written into bench/ but are NOT committed
# (.gitignore excludes them).

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
GS2="${GS2:-$ROOT/gs2}"
INPUT="${INPUT:-$ROOT/example/200.faa}"
E="${E:-100}"
FREQ="${FREQ:-999}"

if [[ ! -x "$GS2" ]]; then
  echo "ERROR: $GS2 not found. Build first: make" >&2
  exit 64
fi
SUB="${1:-stat}"
shift || true
OUT="${OUT:-$ROOT/bench/profile-$SUB}"

have_perf() { command -v perf > /dev/null; }

# Returns 0 if `perf stat` is permitted (paranoid <= 1 or running as root).
perf_allowed() {
  have_perf || return 1
  local p
  p="$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo 4)"
  [[ "$EUID" -eq 0 ]] || [[ "$p" -le 1 ]]
}

case "$SUB" in
  stat)
    if perf_allowed; then
      echo "[profile] perf stat -d  $GS2 -s -e $E -r 12345 $INPUT"
      perf stat -d "$GS2" -s -e "$E" -r 12345 "$INPUT" > /dev/null 2> "$OUT.txt"
      echo "[profile] wrote $OUT.txt"
      cat "$OUT.txt"
    else
      echo "[profile] perf unavailable (paranoid=$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null), euid=$EUID)"
      echo "[profile] falling back to /usr/bin/time -v"
      /usr/bin/time -v -o "$OUT.txt" "$GS2" -s -e "$E" -r 12345 "$INPUT" > /dev/null
      echo "[profile] wrote $OUT.txt"
      cat "$OUT.txt"
      echo
      echo "[profile] to enable perf instead: sudo sysctl kernel.perf_event_paranoid=1"
    fi
    ;;

  record)
    if ! perf_allowed; then
      echo "ERROR: perf record requires perf_event_paranoid<=1 or root." >&2
      echo "  Try: sudo sysctl kernel.perf_event_paranoid=1" >&2
      exit 1
    fi
    echo "[profile] perf record -F $FREQ -g  $GS2 ..."
    perf record -F "$FREQ" -g -o "$OUT.data" -- \
      "$GS2" -s -e "$E" -r 12345 "$INPUT" > /dev/null
    echo "[profile] wrote $OUT.data — run: perf report -i $OUT.data"
    ;;

  flame)
    if ! perf_allowed; then
      echo "ERROR: perf flame requires perf_event_paranoid<=1 or root." >&2
      echo "  Try: sudo sysctl kernel.perf_event_paranoid=1" >&2
      exit 1
    fi
    : "${FLAMEGRAPH_DIR:?Set FLAMEGRAPH_DIR to the path of FlameGraph scripts}"
    perf record -F "$FREQ" -g --call-graph=dwarf -o "$OUT.data" -- \
      "$GS2" -s -e "$E" -r 12345 "$INPUT" > /dev/null
    perf script -i "$OUT.data" \
      | "$FLAMEGRAPH_DIR/stackcollapse-perf.pl" \
      | "$FLAMEGRAPH_DIR/flamegraph.pl" > "$OUT.svg"
    echo "[profile] wrote $OUT.svg"
    ;;

  *)
    echo "usage: $0 {stat|record|flame}" >&2
    exit 64
    ;;
esac
