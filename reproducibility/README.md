# HyperIso reference reproducibility suite

This directory contains the compact reference suite used by the CPC paper.  It is
intended to be distributed with the source archive and to serve as a lightweight
regression check for the public command-line workflows used in the manuscript.

The suite is deliberately small: it does not replace the unit/integration tests.
Instead, it checks that a released executable can reproduce the same Wilson,
observable and statistical summaries from a frozen input file.

## Layout

```text
reproducibility/
  README.md
  manifest.json
  inputs/
    sm_reference.flha
  scripts/
    run_cli_suite.sh
    check_expected_outputs.py
  expected_outputs/
    README.md
  outputs/              # created/overwritten by local runs; not part of the reference archive
```

## Running the suite

Build the command-line executable first, for example

```bash
cmake -S . -B build -DBUILD_WITH_APP=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target hyperiso-ui -j
```

Then run

```bash
./reproducibility/scripts/run_cli_suite.sh
```

The script looks for the executable in a few common build locations and then on
`PATH`.  A custom executable can be forced with

```bash
HYPERISO_BIN=/path/to/hyperiso-ui ./reproducibility/scripts/run_cli_suite.sh
```

The numerical Monte-Carlo examples use `--seed 123456` and `MC_threads=1` by
default.  If the reference outputs are not yet frozen, create them with

```bash
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

Commit the files written in `reproducibility/expected_outputs/` only after the
release metadata, compiler/build configuration and reference input files have
been frozen.

## Comparing with reference outputs

When reference outputs are present, the runner automatically calls
`check_expected_outputs.py`.  The checker compares all floating-point values
extracted from the text outputs with the tolerances declared in `manifest.json`.
The comparison is intentionally tolerant of harmless formatting changes, but it
will fail if a value is missing or outside the declared absolute/relative
threshold.

For Monte-Carlo summaries, the comparison is meaningful only for serial runs
with the fixed seed.  Parallel MC may still be statistically consistent, but it
is not guaranteed to consume the RNG stream in the same order on all platforms.
