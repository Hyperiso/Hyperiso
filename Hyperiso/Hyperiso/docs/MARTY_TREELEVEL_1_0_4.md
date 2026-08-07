# MARTY tree-level matching corrections for HyperIso 1.0.4

This branch corrects the generic MARTY tree-level path used by `Model.MARTY`
and makes the external-fermion ordering configurable per Wilson coefficient and
per perturbative order.

## Per-coefficient external-fermion ordering

MARTY's automatic external-fermion ordering can choose a different pairing for
a four-fermion amplitude at tree level and at one loop. A single global order is
therefore unsafe: it can fix one process while silently breaking another, and it
cannot represent a coefficient whose correct tree and loop projectors differ.

`HyperisoConfig` now exposes two independent maps:

```python
config.mty_tree_fermion_orders = {
    "C9": [1, 0, 2, 3],
    "C_BS_1": [1, 0, 3, 2],
}
config.mty_one_loop_fermion_orders = {
    "C9": [1, 0, 2, 3],
    # This may differ from the tree order when a model/template requires it.
    "C_BS_1": [1, 0, 3, 2],
}
```

An entry affects only the named coefficient and only the selected perturbative
order. Missing entries preserve the order embedded in the corresponding MARTY
template. If a template has no explicit order, MARTY's automatic ordering is
left enabled.

The current four-fermion template defaults are:

- semileptonic and charged-current operators: `[1, 0, 2, 3]`;
- neutral-meson-mixing operators: `[1, 0, 3, 2]`.

The former global `mty_tree_fermion_order` option is removed because it
incorrectly overrode every generated coefficient.

## BSM perturbative-order policy

```python
from pyhyperiso.Common import MartyOrderPolicy

config.mty_order_policy = MartyOrderPolicy.AUTO
```

The available policies are:

- `AUTO`: for a pure BSM coefficient, compute the tree-level contribution first;
  if it is exactly zero, compute the one-loop contribution instead;
- `TREE_LEVEL_ONLY`: retain the tree coefficient even when it is exactly zero;
- `ONE_LOOP_ONLY`: skip the tree probe for a loop-leading validation.

The user-selected policy applies only to the pure BSM target. The independently
generated SM baseline retains its existing automatic behavior. For BSM
coefficients, both templates originally written as TreeLevel-only and templates
originally written as OneLoop-only are routed through the same tree-first
wrapper. Templates that already implement both calls, such as `C10`, keep their
specialised internal fallback. In every BSM one-loop branch, the generated
diagram filter requires at least one non-SM particle, so the result cannot
contain a second copy of the SM matching.

No permutation scan is performed. Exact zeros are physically meaningful, and a
wrong pairing can produce a non-zero but incorrect Fierz projection.

## Generated-code behavior

Both the generic tree-first path and the split-regulator semileptonic path now:

1. resolve the configured order for the current coefficient and perturbative
   order;
2. otherwise retain the explicit order from the coefficient template;
3. disable MARTY's automatic external reordering only when an explicit order is
   present;
4. report the selected perturbative order and fermion order in the generation
   log.

## Kinematic extraction

The generated-source parser recognises wrapped external particles such as
`Outgoing(AntiPart("mu"))` and selects the most complete insertion list instead
of returning a partial 1-to-1 or 1-to-2 process. This removes the legacy KIN
fallback for normal `b -> s l+ l-` templates.

## Neutral-meson 2-to-2 kinematics

`SMParamSetter` now handles the `2 -> 2` external topology used by the
neutral-meson mixing templates.  For elastic mixing processes such as
`b + anti-s -> anti-b + s`, HyperIso evaluates the generated MARTY scalar
products at the on-shell centre-of-mass threshold.  At that point all external
three-momenta vanish and `s_ij = p_i.p_j = m_i m_j`.  This removes the old
fallback to the semileptonic `b -> s mu mu` KIN prescription, which was not
appropriate for Delta F = 2 matching.

The canonical SM mapping now uses `KIN:12` and `KIN:13` for MARTY `s_12` and
`s_13`.  `SMParamSetter` still accepts the historical `KIN:4` and `KIN:7`
aliases so existing user mappings remain compatible.

## ParameterProvider

A typed `ParameterProvider` can be called with a `ParamId` carrying the same type
without emitting a false warning. Untyped identifiers inherit the provider
type; genuine type mismatches are rejected.

## Cache compatibility

The MARTY cache ABI is bumped to `pyhyperiso-1.0.4-v5`. The cache marker includes
the order policy and both the tree-level and one-loop fermion orders for the
coefficient. Existing generated sources and libraries are invalidated
automatically after switching to this branch or changing either map.

## Constructor-time SM inputs

Target MARTY models are constructed only after `undefineNumericalValues()` has
been called. This prevents custom model constructors from folding MARTY's
default values of `theta_W`, `e_em`, or other Standard-Model constants into the
Lagrangian before HyperIso injects the selected LHA inputs.
