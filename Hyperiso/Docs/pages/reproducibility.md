# Reproducibility and release references {#reproducibility}

HyperIso ships a release-level reference suite under `reproducibility/`. The
suite is command-line based so that its inputs, commands and normalized outputs
can be archived and compared without depending on a notebook or graphical
session.

## Reference cases

| ID | Workflow |
|---|---|
| R1 | SM B-sector Wilson coefficients at matching and low scales. |
| R2 | Scalar B-sector Wilson coefficients. |
| R3 | Inclusive and leptonic observables. |
| R4 | A binned angular observable in a fixed \f$q^2\f$ interval. |
| R5 | Fixed-seed Monte-Carlo uncertainty summary and accepted-sample CSV. |
| R6 | THDM Wilson coefficients from an archived 2HDMC spectrum. |
| R7 | SUSY Wilson coefficients from an archived SOFTSUSY spectrum. |

The exact commands, tolerances and output filenames are defined in
`reproducibility/manifest.json`.

## Run the suite

After building `hyperiso-ui` in Release mode:

```bash
./reproducibility/scripts/run_cli_suite.sh
```

Use a non-standard executable with:

```bash
HYPERISO_BIN=/absolute/path/to/hyperiso-ui \
  ./reproducibility/scripts/run_cli_suite.sh
```

The checker normalizes startup banners and absolute paths where they are not
part of the scientific result.

## Statistical determinism

The R5 case uses:

- seed `123456`;
- one Monte-Carlo thread;
- 200 accepted draws;
- a frozen set of output columns and numerical tolerances.

This is a release-to-release regression test, not a convergence claim for a
physics analysis. Production studies should choose a draw count based on their
own precision requirements.

## Frozen provenance

Reference outputs are regenerated only for an intentional numerical change:

```bash
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

The freeze records:

- the source commit used to build the executable;
- the executable SHA-256 hash;
- compiler, CMake, GSL, Eigen and Python versions;
- build type and platform information;
- per-case timings;
- hashes of the expected outputs.

Review every reference diff together with the physics or numerical change that
caused it. The provenance commit should contain only the frozen expected-output
files.

## Reproducing external-model cases

The THDM and SUSY release checks consume archived spectrum files. They do not
rerun 2HDMC or SOFTSUSY. This keeps the mandatory release gate independent of
external-generator availability while retaining the numerical model input.

MARTY validation is documented separately and is not part of the mandatory
frozen suite because it requires an external symbolic installation.

## Release identity

For a published result, retain:

- the immutable Git tag and commit;
- the source archive and checksums;
- the exact input and override files;
- the software environment or wheel filename;
- the reference manifest and provenance metadata;
- the random seed and thread configuration.
