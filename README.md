# PP: Pairwise Sequence Alignment Benchmark

This project benchmarks and compares two C-based pairwise sequence alignment implementations:

- a rolling-row Smith-Waterman implementation in `src/tiled_sw.c`
- a banded dynamic-programming implementation in `src/wfa_banded.c`

The project also includes a Python harness for generating synthetic sequence pairs, running the binaries across multiple thread counts and identity/length settings, and collecting timing and memory statistics.

## Project structure

- `src/` — C alignment implementations
- `harness/` — benchmark orchestration script
- `results/` — generated FASTA inputs, CSV summaries, and plots

## Benchmark workflow

1. Compile the C binaries.
2. Generate benchmark inputs with the Python harness.
3. Each invocation runs one requested thread count and dynamic-schedule chunk.
4. Inspect raw/aggregated results and scaling plots in `results/`.

The binaries can also be run directly:

```text
tiled_sw.exe pairs.fa [threads [schedule_chunk]]
wfa_banded.exe pairs.fa band [threads [schedule_chunk]]
```

## Example build

From the repository root:

```bash
gcc -O3 -fopenmp -o tiled_sw.exe src/tiled_sw.c
gcc -O3 -fopenmp -o wfa_banded.exe src/wfa_banded.c
```

With `--auto-compile`, the harness tries GCC, MSVC (`cl`), or Clang and uses
the matching OpenMP compiler flags. Select one explicitly with
`--compiler gcc`, `--compiler clang`, or `--compiler cl` when needed. MSVC
must be launched from a Visual Studio Developer Command Prompt with OpenMP
support installed.

## Example harness run

```bash
python harness/run_experiments_harness.py \
  --outdir results \
  --pairs-per-identity 200 \
  --lengths 2000 5000 \
  --identities 0.98 0.95 0.90 0.80 \
  --threads-list 1 2 4 8 \
  --schedule-chunks 1 2 4 \
  --repeats 3 \
  --auto-compile
```

## Linux hardware profiling

On Linux, `--perf` wraps each run in `perf stat` and records L1 data-cache
load misses, instructions, cycles, and derived IPC in the CSV files. Event
availability and perf permissions depend on the host CPU and kernel.

CPU affinity is optional and uses Linux `taskset`. Supply CPU lists for the
topology on the machine being tested; CPU IDs are not portable between
machines:

```bash
python harness/run_experiments_harness.py \
  --threads-list 1 2 4 8 \
  --schedule-chunks 1 2 4 \
  --repeats 3 \
  --perf \
  --p-core-cpus 0-5 \
  --e-core-cpus 6-11 \
  --all-core-cpus 0-11
```

The raw and aggregated CSVs include the affinity profile and CPU set, and the
speedup/scaling plots separate the profiles and chunk sizes. The harness
reports measured results; it does not assume or guarantee linear scaling.

## Windows affinity and Intel VTune

Windows runs can restrict each benchmark process to logical processors using
CPU affinity masks. Determine the P-core and E-core logical processor IDs on
the machine first, then create a mask with one bit set per selected processor.
For example, `0x3` selects logical processors 0 and 1; these example masks
are not a claim about any particular Core Ultra topology.

```powershell
python .\harness\run_experiments_harness.py `
  --threads-list 1 2 4 8 `
  --schedule-chunks 1 2 4 `
  --p-core-mask 0x3 `
  --e-core-mask 0xC `
  --all-core-mask 0xF
```

Windows masks are applied directly to the benchmark process with
`SetProcessAffinityMask` and are inherited by its OpenMP workers. The harness
applies the mask immediately after process creation. This mode is limited to
the process's Windows processor group; systems spanning multiple processor
groups need a group-aware affinity implementation.

To collect Intel VTune microarchitecture data for each run, add `--vtune`.
The default analysis is `uarch-exploration`; it can be changed with
`--vtune-analysis`. If VTune's `bin64` directory is not already on `PATH`,
pass it using `--vtune-bin`. That directory is prepended only to the harness's
child-process `PATH`; the machine-wide System PATH is not modified.

```powershell
python .\harness\run_experiments_harness.py `
  --vtune `
  --vtune-bin "C:\Program Files (x86)\Intel\oneAPI\vtune\latest\bin64" `
  --p-core-mask 0x3 `
  --e-core-mask 0xC `
  --all-core-mask 0xF
```

The harness exports a VTune CSV summary per run and parses recognized L1 data
cache-miss, retired-instruction, cycle, and IPC metrics into the raw and
aggregated CSVs. Unrecognized or unavailable report fields remain empty; the
VTune result and exported summary paths are retained for manual inspection.
Metric names and collection support vary with VTune version, CPU, and
permissions. Do not enable `--perf` and `--vtune` together.

## Notes

- `tiled_sw.c` is currently a rolling-row implementation, not a tiled
  wavefront implementation. Both algorithms parallelize over independent
  sequence pairs using OpenMP dynamic scheduling.
- The harness expects the compiled binaries to be available in the project
  root or in the paths supplied to it.
- Results are stored in CSV and plotting output under `results/`.
- The Linux `perf` and `taskset` options are opt-in; default benchmark runs do
  not require those tools.

## License

This project is distributed under the MIT license. See the `LICENSE` file for details.
