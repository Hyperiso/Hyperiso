# Frozen HyperIso 1.0.0 reference outputs

These files are normalized CPC reference results. The startup banner and all
machine-specific paths are intentionally excluded. `statistics_samples.csv`
contains exactly 200 accepted serial Monte-Carlo samples generated with seed
`123456`.

Update the references only through:

```bash
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

The command validates finite values, CSV structure and sample count before
replacing the references. Every update must be reviewed as a numerical change.
All CSV fields are frozen with the full precision emitted by the 1.0.0 writer. Text and CSV comparisons account only for the decimal resolution represented by each frozen token, so a last-digit formatting boundary is not mistaken for a physics change.
