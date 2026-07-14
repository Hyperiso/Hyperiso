# HyperIso reference reproducibility suite

This directory contains the compact reference suite used for the CPC release.
It checks the public CLI workflows against frozen, machine-independent outputs.
It complements the unit and integration tests; it does not replace them.

## Layout

```text
reproducibility/
  manifest.json
  inputs/sm_reference.flha
  scripts/
    run_cli_suite.sh
    normalize_cli_output.py
    check_expected_outputs.py
    freeze_reference_outputs.py
  expected_outputs/
  outputs/                  # local, ignored working directory
```

## Run

Build the release CLI first:

```bash
cmake -S Hyperiso/Hyperiso/core -B build \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_WITH_CLI=ON
cmake --build build --parallel
```

Then run:

```bash
./reproducibility/scripts/run_cli_suite.sh
```

A custom executable can be selected with:

```bash
HYPERISO_BIN=/path/to/hyperiso-ui \
  ./reproducibility/scripts/run_cli_suite.sh
```

The statistical case uses **1,000 accepted serial draws** and seed `123456`.
The generated text files begin at the stable `== ... summary ==` section; the
startup banner and machine-specific paths are deliberately excluded.

## Updating references

Only update references after reviewing a deliberate numerical change:

```bash
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

The update is rejected when a text output is missing or non-normalized, when any
value is non-finite, or when the CSV does not contain exactly 1,000 finite rows
with the declared columns. The metadata file records SHA-256 hashes for every
frozen reference artifact.

The comparison uses both absolute and relative tolerances through
`math.isclose`; missing references are fatal. Parallel Monte-Carlo execution is
not used for reference generation because worker scheduling can change random
number consumption order.
