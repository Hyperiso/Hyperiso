# Frozen HyperIso 1.0.0 reference outputs

These files are normalized CPC reference results. The startup banner and all
machine-specific paths are intentionally excluded. `statistics_samples.csv`
contains exactly 1,000 accepted serial Monte-Carlo samples generated with seed
`123456`.

Update the references only through:

```bash
./reproducibility/scripts/run_cli_suite.sh --update-expected
```

The command validates finite values, CSV structure and sample count before
replacing the references. Every update must be reviewed as a numerical change.
