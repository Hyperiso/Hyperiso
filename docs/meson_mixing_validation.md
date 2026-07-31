# Neutral-meson-mixing validation and conventions

The `MesonMixing` implementation has been validated through independent tests of
matching, QCD running, basis transformations, hadronic assembly and observables.
The static guard is:

```bash
python tools/check_meson_mixing_regressions.py
```

## Validated chain

For the public SUSY-basis coefficients

```text
C1, C2, C3, C4, C5, CT1, CT2, CT3
```

in the `BD`, `BS`, `SD` and `CU` families, the following layers were checked:

1. real and complex `FWCOEF`/`IMFWCOEF` input;
2. five-flavour running for `BD/BS` and four-flavour running for `SD/CU`;
3. NLO expansion with `alpha_s(mu_low)/(4*pi)`;
4. SUSY -> BMU -> SUSY transformations;
5. family-specific scales `B_SCALE`, `K_SCALE`, `D_SCALE`;
6. matrix elements `Q1`--`Q5`, including the chiral factor for `Q2`;
7. `M12` assembly and the exposed observables;
8. tree-level neutral-vector matching:

   ```text
   C1  = Delta_L^2 / (2 M^2)
   CT1 = Delta_R^2 / (2 M^2)
   C5  = -2 Delta_L Delta_R / M^2
   ```

The direct MARTY chain and an independent analytic `FWCOEF/IMFWCOEF` chain agree
at the `1.2e-6` level in matching, running and final observables for the generic
mixed-chiral Z-prime benchmark at 3 TeV.


## Separation of SM and BSM matching scales

`EW_SCALE` is the matching scale of the Wilson-coefficient group and may be set
to a heavy mediator mass.  The Standard-Model W-box contribution in
`M0_Mixing` must not use that BSM scale.  Its electroweak scale is instead tied
to `m_W` and the existing `SCALE_NUIS[1]` nuisance.

Consequently, varying a BSM matching scale from the electroweak scale to a
multi-TeV scale must leave all SM-only meson-mixing observables invariant.

## Observable conventions and limitations

* `PHI_D` is the historical name of the **Bd mixing-amplitude phase**, not the
  D-meson phase. Do not compare it directly with a D-mixing phase.
* `PHI_S` is the phase of the Bs mixing amplitude `arg(M12)`. It is not by
  itself the decay-interference phase `phi_s^{ccs}` measured in `Bs -> J/psi phi`.
* `A_FS_D` and `A_FS_S` assume new physics modifies `M12` but not `Gamma12`.
* `DELTA_M_D` and `X_D` currently contain the short-distance NP contribution
  only; the long-distance SM D-mixing amplitude is not modelled.
* `DELTA_M_K` contains the implemented short-distance SM term plus NP, not a
  complete long-distance SM prediction. It must not be used as a naive precise
  Gaussian constraint on the experimental total.
* `ABS_EPSILON_K` uses the experimental `Delta M_K` denominator and is the more
  suitable kaon observable for a first phenomenological scan.

## Parameter numbering

The current `M0_Mix` convention is:

```text
5_1,5_2 : |Gamma12_s|, arg(Gamma12_s)
6       : DeltaGamma_s
7_1,7_2 : |Gamma12_d|, arg(Gamma12_d)
8       : DeltaGamma_d
9       : kappa_epsilon
10      : mu_B
11      : mu_K
12      : mu_D
13      : eta_cc
14      : eta_ct
15      : experimental DeltaM_K in s^-1
```
