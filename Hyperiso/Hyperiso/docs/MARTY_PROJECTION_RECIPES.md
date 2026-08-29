# MARTY projection recipes

## Why recipes exist

MARTY four-fermion matching has several independent choices:

- `F`: external/Fierz fermion order used to construct and match the amplitude;
- `O`: fermion pairing used by `dimension6Operator`;
- Dirac-current basis (`VL x V`, `VL x A`, `VL x VL`, ...);
- current layout (quark current first or lepton current first);
- linear combinations of several primitive projections.

No universal rule maps chirality alone to all of these choices.  In particular,
the validated direct neutral-vector Z' projection and the validated t-channel LQ
projection use different recipes.

Normal HyperIso calculations **never scan** these choices.  They either use the
coefficient template/default F/O settings or a frozen explicit recipe.

## Precedence

For TreeLevel `C9`, `C10`, `CP9`, and `CP10`:

1. an explicit projection recipe owns F/O/current/layout term-by-term;
2. otherwise the existing `mty_tree_fermion_orders` and
   `mty_tree_operator_orders` entries are used;
3. otherwise the historical template defaults are used.

OneLoop behavior is unchanged by recipe ABI v1.

All other Wilson coefficients keep the existing generic per-coefficient F/O
configuration unchanged.

## Python API

```python
from pyhyperiso.marty import (
    MartyProjectionProfile,
    apply_tree_projection_profile,
    direct_semileptonic_profile,
    lq_semileptonic_chiral_profile,
)
```

### Direct neutral-vector projection

The default template is direct, so a Z' workspace may simply keep its existing
per-coefficient F/O maps.  An explicit frozen profile can also be used:

```python
profile = direct_semileptonic_profile(
    fermion_order=[0, 2, 1, 3],
    operator_order=[0, 2, 1, 3],
)
apply_tree_projection_profile(config, profile)
```

This encodes

```text
C9   : VL x V
C10  : VL x A
CP9  : VR x V
CP10 : VR x A
```

### Validated LQ chiral projection

```python
profile = lq_semileptonic_chiral_profile()
apply_tree_projection_profile(config, profile)
```

This freezes

```text
LL : F=[1,0,2,3], O=[0,3,1,2]
LR : F=[1,0,2,3], O=[1,2,0,3]
RL : F=[1,0,2,3], O=[1,2,0,3]
RR : F=[1,0,2,3], O=[1,2,0,3]
```

and reconstructs

```text
C9   = +0.5 LL +0.5 LR
C10  = -0.5 LL +0.5 LR
CP9  = +0.5 RL +0.5 RR
CP10 = -0.5 RL +0.5 RR
```

The helper removes only coarse TreeLevel C9/C10/CP9/CP10 F/O entries that would
otherwise conflict with the term-by-term recipe. OneLoop entries and unrelated
coefficients are preserved.

## Profile files

Profiles can be frozen as JSON:

```python
profile.save("marty_projection_profile.json")
profile = MartyProjectionProfile.load("marty_projection_profile.json")
```

This is the recommended form for paper/reproducibility workspaces.

## Diagnostic scanner

`SubprocessProjectionScanner` evaluates candidates in isolated subprocesses.
It does not modify the production HyperIso process and it never silently selects
the first non-zero permutation.

A workspace runner needs only two hooks:

```python
from pyhyperiso.marty import (
    apply_projection_profile_from_env,
    write_projection_scan_result,
)

config = make_my_hyperiso_config()
apply_projection_profile_from_env(config)

# initialize HyperIso and calculate the requested matching coefficients
values = {"C9": complex(c9)}
write_projection_scan_result(values)
```

Then a diagnostic script can run e.g.

```python
from pyhyperiso.marty import SubprocessProjectionScanner, all_fermion_orders

scanner = SubprocessProjectionScanner(
    ["python", "checks/projection_runner.py", "--model", "my_model"],
    cwd="/path/to/workspace",
)

results = scanner.scan_single_term(
    coefficient="C9",
    fermion_orders=all_fermion_orders(),
    operator_orders=all_fermion_orders(),
    current_pairs=[("VL", "V"), ("VL", "A"), ("VL", "VL"), ("VL", "VR")],
    layouts=["quark_first", "lepton_first"],
    expected=0.4512083979756 + 0j,  # optional analytic/reference value
)

for candidate in results[:10]:
    print(candidate.score, candidate.value, candidate.profile.to_json())
```

An unrestricted 24 x 24 x currents x layouts scan can be expensive.  Restrict
candidate families whenever physics or a standalone MARTY amplitude diagnostic
already determines part of the layout.

The scanner sorts by relative error only when an explicit reference value is
provided.  Without a reference it reports candidates by magnitude and leaves
selection to the validation workflow.  It never treats "first non-zero" as a
physical criterion.

## Cache

Recipe ABI v1 is part of the generated-source cache marker.  Changing any term's
weight, F, O, currents or layout invalidates the affected generated coefficient.
The MARTY cache ABI is bumped to `pyhyperiso-1.0.4-v12`.

### Command-line scan

The same scanner is exposed as a CLI:

```bash
python -m pyhyperiso.marty.scan_projection \
  --runner "python checks/projection_runner.py --model my_model" \
  --coefficient C9 \
  --fermion-order 0,2,1,3 \
  --fermion-order 1,0,2,3 \
  --operator-order 0,2,1,3 \
  --operator-order 0,3,1,2 \
  --current-pair VL,V \
  --current-pair VL,VL \
  --current-pair VL,VR \
  --all-layouts \
  --expected-real 0.4512083979756 \
  --top 20 \
  --export-best marty_projection_profile.json
```

Use `--all-fermion-orders` and/or `--all-operator-orders` for exhaustive
permutation scans.  Candidate subprocess failures are recorded and skipped so a
single invalid MARTY pairing does not abort the whole diagnostic.  Exporting the
largest non-zero result without a reference is refused by default; that safety
check can only be bypassed explicitly with `--allow-magnitude-ranking`.

For a new model, the recommended validation sequence is:

1. inspect the standalone MARTY amplitude and restrict obvious impossible
   layouts;
2. scan candidate F/O/current/layout combinations;
3. compare against an analytic value or independent implementation when one is
   available;
4. vary at least one coupling and one mediator mass to test scaling/stability;
5. test linearity with mixed textures when several primitive currents can be
   active simultaneously;
6. freeze the accepted JSON profile and use only that profile in scans/fits.
