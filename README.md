# PP: Pairwise Sequence Alignment Benchmark

This project benchmarks and compares two C-based pairwise sequence alignment implementations:

- a tiled Smith-Waterman implementation in `src/tiled_sw.c`
- a banded WFA-style implementation in `src/wfa_banded.c`

The project also includes a Python harness for generating synthetic sequence pairs, running the binaries across multiple thread counts and identity/length settings, and collecting timing and memory statistics.

## Project structure

- `src/` — C alignment implementations
- `harness/` — benchmark orchestration script
- `results/` — generated FASTA inputs, CSV summaries, and plots

## Benchmark workflow

1. Compile the C binaries.
2. Generate benchmark inputs with the Python harness.
3. Run the alignment binaries across the desired parameter sweep.
4. Inspect the raw and aggregated results in `results/`.

## Example build

From the repository root:

```bash
gcc -O3 -fopenmp -o tiled_sw src/tiled_sw.c
gcc -O3 -fopenmp -o wfa_banded src/wfa_banded.c
```

## Example harness run

```bash
python harness/run_experiments_harness.py \
  --outdir results \
  --pairs-per-identity 200 \
  --lengths 2000 5000 \
  --identities 0.98 0.95 0.90 0.80 \
  --threads-list 1 2 4 8 \
  --repeats 3 \
  --auto-compile
```

## Notes

- The harness expects the compiled binaries to be available in the project root or in the paths supplied to it.
- Results are stored in CSV and plotting output under `results/`.
- The implementation focuses on benchmark automation and performance comparison rather than a production-ready CLI.

## License

This project is distributed under the MIT license. See the `LICENSE` file for details.
