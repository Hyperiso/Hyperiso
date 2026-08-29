from types import SimpleNamespace

from pyhyperiso.marty.projection import (
    TREE_RECIPE_PREFIX,
    MartyProjectionProfile,
    apply_tree_projection_profile,
    direct_semileptonic_profile,
    lq_semileptonic_chiral_profile,
)


def test_lq_profile_encodes_four_coefficients_and_preserves_unrelated_orders():
    cfg = SimpleNamespace(
        mty_tree_fermion_orders={"C7": [1, 0, 2, 3], "C9": [3, 2, 1, 0]},
        mty_tree_operator_orders={"C7": [1, 0, 2, 3], "C9": [3, 2, 1, 0]},
    )
    apply_tree_projection_profile(cfg, lq_semileptonic_chiral_profile())
    assert cfg.mty_tree_fermion_orders["C7"] == [1, 0, 2, 3]
    assert cfg.mty_tree_operator_orders["C7"] == [1, 0, 2, 3]
    assert "C9" not in cfg.mty_tree_fermion_orders
    assert "C9" not in cfg.mty_tree_operator_orders
    keys = [key for key in cfg.mty_tree_fermion_orders if key.startswith(TREE_RECIPE_PREFIX)]
    assert len(keys) == 8
    assert all(key in cfg.mty_tree_operator_orders for key in keys)


def test_direct_profile_roundtrip(tmp_path):
    profile = direct_semileptonic_profile([0, 2, 1, 3], [0, 2, 1, 3])
    path = profile.save(tmp_path / "profile.json")
    loaded = MartyProjectionProfile.load(path)
    assert loaded.to_json() == profile.to_json()
