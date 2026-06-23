# Reproducibility plan

This page describes the recommended structure for reproducing the numerical results shown in the HyperIso paper.

## Goals

The reproducibility package should allow users to regenerate:

- Wilson-coefficient comparison tables;
- observable validation plots;
- SM benchmark predictions;
- MARTY-vs-native model checks;
- statistical uncertainty and fit examples.

## Proposed structure

```text
reproducibility/
├── README.md
├── inputs/
│   ├── sm/
│   ├── thdm/
│   ├── susy/
│   └── marty/
├── scripts/
│   ├── run_cpp_examples.sh
│   ├── run_python_examples.py
│   ├── run_cli_examples.sh
│   └── make_validation_plots.py
├── expected/
│   ├── wilson_sm.json
│   ├── wilson_thdm.json
│   ├── wilson_susy.json
│   └── observables_sm.json
└── figures/
```

## Numerical-output policy

For each benchmark, store:

- input file path;
- command;
- HyperIso commit hash;
- backend flags;
- numerical result;
- absolute and relative tolerance;
- reference source.

## CI policy

Fast reproducibility tests should run on every pull request. Heavy validation should run on schedule or through `workflow_dispatch`.
