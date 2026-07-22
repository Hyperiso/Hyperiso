# MARTY tree-level matching corrections for HyperIso 1.0.4

This branch corrects the generic MARTY tree-level path used by `Model.MARTY`.

## Correct external-fermion ordering

MARTY's three-argument `computeWilsonCoefficients(...)` overload enables its
automatic external-fermion ordering. In a four-fermion process, that automatic
pairing may differ between tree and one-loop amplitudes and may override an
explicit `FeynOptions::setFermionOrder(...)` supplied by a coefficient template.

HyperIso 1.0.4 now preserves an explicit template order for every tree-level
Wilson calculation. Templates that do not provide a fermion order retain
MARTY's automatic ordering, and one-loop calculations keep the historical
behavior. This applies to ordinary coefficient templates, the generic
tree-first wrapper, and the split-regulator semileptonic path.

## Explicit order policy

`HyperisoConfig` now exposes:

```python
from pyhyperiso.Common import MartyOrderPolicy

config.mty_order_policy = MartyOrderPolicy.TREE_LEVEL_ONLY
```

The available policies are:

- `AUTO`: historical tree-first, then one-loop fallback when the tree
  coefficient vanishes;
- `TREE_LEVEL_ONLY`: retain the tree coefficient even when it is exactly zero;
- `ONE_LOOP_ONLY`: skip the tree probe for a loop-leading validation.

The policy is applied only to the configured BSM target. The separately
generated Standard-Model baseline remains `AUTO`.

## Kinematic extraction

The generated-source parser now recognises wrapped external particles such as
`Outgoing(AntiPart("mu"))` and selects the most complete insertion list instead
of returning a partial 1-to-1 or 1-to-2 process. This removes the legacy KIN
fallback for normal `b -> s l+ l-` templates.

## ParameterProvider

A typed `ParameterProvider` can now be called with a `ParamId` carrying the same
type without emitting a false warning. Untyped identifiers inherit the
provider type; genuine type mismatches are rejected.

## Cache compatibility

The MARTY cache ABI is bumped to `pyhyperiso-1.0.4-v1`, and the configured order
policy is part of the generation-mode marker. Existing generated sources and
libraries are therefore invalidated automatically after switching to this
branch.

## Why HyperIso does not scan permutations automatically

The fix preserves an order explicitly requested by a coefficient template; it
does not select the first non-zero result among all permutations. Exact zeros
are physically meaningful, and a wrong pairing can also give a non-zero but
incorrect Fierz projection. The coefficient/model implementation must therefore
provide the intended fermion order. If it does not, MARTY's normal automatic
ordering remains enabled.
