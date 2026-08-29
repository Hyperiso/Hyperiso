# MARTY TreeLevel semileptonic chiral matching

> **Superseded implementation note (cache ABI v11).**  The LQ chiral
> decomposition is no longer hard-coded as the universal TreeLevel path for
> `C9`, `C10`, `CP9`, and `CP10`.  It is retained as the built-in
> `lq_semileptonic_chiral_profile()` projection recipe.  See
> `MARTY_PROJECTION_RECIPES.md` for the general mechanism.

The validated LQ decomposition remains

    C9   = (LL + LR) / 2
    C10  = (-LL + LR) / 2
    CP9  = (RL + RR) / 2
    CP10 = (-RL + RR) / 2

with

    F(LL/LR/RL/RR) = [1, 0, 2, 3]
    O(LL)           = [0, 3, 1, 2]
    O(LR/RL/RR)     = [1, 2, 0, 3]

These values are a validated **profile**, not a universal MARTY theorem.  A
different topology or external-current layout may require another recipe.  In
particular, direct neutral-vector exchange can use the direct `VL x V/A` or
`VR x V/A` projectors instead of the LQ chiral decomposition.

The generic per-coefficient controls remain available for every Wilson:

- `mty_tree_fermion_orders`;
- `mty_one_loop_fermion_orders`;
- `mty_tree_operator_orders`;
- `mty_one_loop_operator_orders`.

An explicit projection recipe is authoritative for the TreeLevel coefficient it
covers and supplies its own F/O term by term.  All OneLoop controls are
unchanged.  HyperIso does not scan permutations during production calculations;
projection discovery is an explicit Python diagnostic step and the accepted
profile should be frozen for scans and fits.
