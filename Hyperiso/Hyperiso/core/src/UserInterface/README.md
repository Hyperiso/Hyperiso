# Hyperiso terminal interface

This folder contains the lightweight terminal interface for common Hyperiso summaries.
It intentionally focuses on stable builtin workflows: Wilson summaries, observable
summaries and Statistic dependency/uncertainty summaries.  Custom runtime decays and
custom observables remain examples/API features rather than CLI primitives.

## Build

By default the library is built but the executable is disabled, following the project-wide
`BUILD_WITH_APP` convention.

```bash
cmake -S . -B build -DBUILD_WITH_APP=ON
cmake --build build --target hyperiso-ui -j
```

From an installed build, replace `./build/hyperiso-ui` by `hyperiso-ui` in the commands below.

## Common options

```bash
--model SM|THDM|MSSM|MARTY      # default: SM
--lha lha/si_input.flha         # default: lha/si_input.flha
--order LO|NLO|NNLO             # default: NNLO
```

## Wilson commands

### Default B-coefficient summary

```bash
./build/hyperiso-ui wilson summary \
  --model SM \
  --lha lha/si_input.flha
```

### Specific Wilson coefficients

```bash
./build/hyperiso-ui wilson summary \
  --groups BCoefficients \
  --coeffs C7,C9,C10 \
  --qmatch 81 \
  --q 4.8 \
  --order NNLO
```

### Several builtin groups

```bash
./build/hyperiso-ui wilson summary \
  --groups BCoefficients,BPrimeCoefficients,BScalarCoefficients \
  --coeffs C7,C9,C10 \
  --order LO
```

## Observable commands

### Default observable summary

```bash
./build/hyperiso-ui observable summary \
  --model SM \
  --lha lha/si_input.flha
```

### Select observables by canonical mapper name

```bash
./build/hyperiso-ui observable summary \
  --observables BR_Bs__mu_mu,BR_B__Xs_gamma \
  --order NNLO
```

### Include a binned observable

`--bins` uses the format `OBS:min:max`.  Multiple bins can be comma-separated.

```bash
./build/hyperiso-ui observable summary \
  --observables BR_Bs__mu_mu \
  --bins F_L_B0__K*0_mu_mu:1.1:6.0 \
  --order NNLO
```

## Statistic commands

### Prediction and dependency summary

```bash
./build/hyperiso-ui statistic summary \
  --observables BR_Bs__mu_mu,BR_B__Xs_gamma \
  --draws 200
```

### Use the chi-square covariance backend

```bash
./build/hyperiso-ui statistic summary \
  --observables BR_Bs__mu_mu,BR_B__Xs_gamma \
  --chi2 \
  --draws 500
```

### Compute MC uncertainty summaries with a progress bar

```bash
./build/hyperiso-ui statistic summary \
  --observables BR_Bs__mu_mu,BR_B__Xs_gamma \
  --uncertainties \
  --progress \
  --draws 500
```

### Verbose diagnostic summary

```bash
./build/hyperiso-ui statistic summary \
  --observables BR_Bs__mu_mu \
  --verbose \
  --draws 100
```

## Help

```bash
./build/hyperiso-ui --help
./build/hyperiso-ui wilson --help
./build/hyperiso-ui observable --help
./build/hyperiso-ui statistic --help
```
