# bench/ — gs2 benchmarking & regression harness

This directory holds scripts for measuring `gs2` wall-time and for guarding
against algorithmic regressions during refactors. It is not built by the
top-level Makefile; the scripts assume `./gs2` has been built (`make`).

The committed baseline corresponds to the `7b317fe` ("v2.5") build with the
bundled reference BLAS. Re-generate it (`bench/baseline/200_e100_seed12345.nwk`)
whenever a PR intentionally changes the algorithm or the BLAS implementation.

## ⚠️ What the bench measures

`run_bench.sh` times the **whole `gs2` invocation, including MMseqs2 / BLASTN
alignment**. For the bundled `example/200.faa` (200 protein sequences),
MMseqs2 dominates by a large margin — on a 64-core EPYC 9554, the four runs
`-e {0, 10, 50, 100}` all land within a few hundredths of ~11.5 s wall, meaning
GS + EP themselves are below ~0.5 s.

That makes the default bench useful for **end-to-end regression checking** but
**near-useless for measuring speedups in `spectral_clustering` / `dsyevr` / `EP_*`**.
To isolate gs2 internals, run gs2 once, capture an `_mmseqs.txt` (or pre-build
a `-w`-format similarity matrix), then time the pre-aligned path — a small
helper for this will land alongside PR-2 where the gs2-internal optimisations
start. Concurrent gs2 invocations on the same input also collide on the shared
`example/200_tmp/` directory (gs2 derives intermediate paths from the input
filename); avoid running bench scripts in parallel against the same FASTA.

## Files

| File | What it does |
|---|---|
| `baseline/200_e100_seed12345.nwk` | Reference Newick for `gs2 -s -e 100 -r 12345 example/200.faa` produced under `OMP_NUM_THREADS=1`. Bit-deterministic. |
| `regression.sh` | Runs gs2 against the baseline. `--strict` = exact diff (OMP=1). Default = topology must match + EP support within `--eps` (0.10). |
| `compare_newick.py` | Pure-stdlib Newick parser; compares two trees by clade set + per-clade EP support tolerance. |
| `run_bench.sh` | Wall-clock benchmark across several `-e` values; reports min/median/max over N repeats. |
| `profile.sh` | Wraps `perf stat` / `perf record` (`stat`/`record`/`flame` subcommands). |

## Determinism

Under `OMP_NUM_THREADS=1`, `gs2 -r 12345` is **bit-deterministic** across runs.
With multiple OpenMP threads, topology stays identical but EP support values
jitter by ≲ 0.06 (typical at 100 iters / 200 seqs) because per-thread `local_ep`
maps merge into the shared `ep` map in scheduler-dependent order, and FP
accumulation isn't associative.

That's why we ship two regression modes:

- **`--strict`** (default) — bit-exact diff under `OMP_NUM_THREADS=1`. Use this
  for PRs that should preserve FP order (`z[N*N]` → `z[N]`, scratch buffer
  reuse, etc.).
- **`--tolerant`** — clade-set equality + `|EP_a − EP_b| ≤ eps` per clade. Use
  this for PRs that legitimately reorder accumulation or otherwise produce
  statistically-equivalent EP draws (alternative eigensolver, OMP scheduling
  restructure, RNG engine swap).

Choosing `--eps`: at `ep_num=100` the binomial noise on each clade's EP is
roughly `1/√100 = 0.10`, with occasional 2σ outliers up to ~0.20. Use **0.20**
for routine tolerant checks, and tighten to **0.10** only when bumping `ep_num`
to ≥ 400.

## Typical use

```bash
# Before changing anything: sanity-check the baseline reproduces.
bench/regression.sh

# After a refactor:
make && bench/regression.sh                            # FP-preserving change
make && bench/regression.sh --tolerant --eps 0.20      # numerical change

# Quick wall-time check after a perf change:
bench/run_bench.sh
THREADS=1 bench/run_bench.sh
REPEATS=5 E_VALUES="0 100 500" bench/run_bench.sh

# Profiling (output files are not committed; .gitignore excludes them):
bench/profile.sh stat                # quick perf stat summary
bench/profile.sh record              # then run: perf report -i bench/profile-record.data
FLAMEGRAPH_DIR=~/FlameGraph bench/profile.sh flame
```

## Regenerating the baseline

```bash
OMP_NUM_THREADS=1 ./gs2 -s -e 100 -r 12345 example/200.faa \
  > bench/baseline/200_e100_seed12345.nwk
```

Commit the regenerated file in the same PR that justified the change, and call
out in the PR description that the baseline was rotated (so reviewers know to
re-run `--strict` against the *new* baseline, not the old one).
