# Reproducibility

HyperIso ships a compact CPC reference reproducibility suite in the repository
root under `reproducibility/`.  The suite is the canonical place for the input
files, commands, frozen reference outputs and numerical tolerances used to check
that the public command-line workflows reported in the paper can be regenerated
from the released source archive.

Run it with

```bash
./reproducibility/scripts/run_cli_suite.sh
```

or, for a non-standard executable location,

```bash
HYPERISO_BIN=/path/to/hyperiso-ui ./reproducibility/scripts/run_cli_suite.sh
```

Reference outputs should only be updated when the release tag, input files,
compiler/build setup and numerical tolerances are frozen:

```bash
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

The Monte-Carlo statistic example uses a fixed `--seed` and the serial MC path by
default.  Parallel MC can be used for performance studies, but serial seeded runs
are recommended for release-to-release numerical comparison.
