# Reproducibility

HyperIso ships a fixed-seed CPC reference suite under `reproducibility/`.
It contains one frozen FLHA input, five normalized CLI summaries and a
200-row Monte-Carlo sample table.

Run the suite after building `hyperiso-ui` in `Release` mode:

```bash
./reproducibility/scripts/run_cli_suite.sh
```

Use a non-standard executable with:

```bash
HYPERISO_BIN=/path/to/hyperiso-ui \
  ./reproducibility/scripts/run_cli_suite.sh
```

Reference outputs intentionally exclude the startup banner and absolute paths,
so comparisons do not depend on the checkout location or user cache directory.
The statistical case uses seed `123456`, one thread and 200 accepted draws.

References may be regenerated only for an intentional numerical update:

```bash
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

The command validates the outputs before freezing them, then records SHA-256
hashes in `expected_outputs/reference_metadata.json`. Review every reference
diff together with the corresponding physics or numerical change.
