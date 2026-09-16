from pathlib import Path
from types import SimpleNamespace

from pyhyperiso.marty.projection import (
    TREE_RECIPE_PREFIX,
    MartyProjectionProfile,
    MartyProjectionRecipe,
    MartyProjectionTerm,
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


def test_bnunu_profile_accepts_direct_neutrino_projector():
    profile = MartyProjectionProfile(
        tree={
            "CNU_L_MUMU": MartyProjectionRecipe(
                "CNU_L_MUMU",
                (MartyProjectionTerm(
                    "direct", 1.0, "VL", "VL",
                    [0, 1, 2, 3], [1, 2, 0, 3], "quark_first"
                ),),
            )
        }
    )
    cfg = SimpleNamespace(mty_tree_fermion_orders={}, mty_tree_operator_orders={})
    apply_tree_projection_profile(cfg, profile)
    keys = [key for key in cfg.mty_tree_fermion_orders if "CNU_L_MUMU" in key]
    assert len(keys) == 1
    key = keys[0]
    assert cfg.mty_tree_fermion_orders[key] == [0, 1, 2, 3]
    assert cfg.mty_tree_operator_orders[key] == [1, 2, 0, 3]


def test_semileptonic_recipe_templates_reorder_only_nontrivial_fermion_orders():
    template_dir = Path(__file__).resolve().parents[2] / "assets" / "template" / "MARTY"
    for name in ("C9.cpp", "C10.cpp", "CP9.cpp", "CP10.cpp"):
        source = (template_dir / name).read_text()
        start = source.index("Expr hyperiso_marty_project_tree_recipe(")
        end = source.index("} // namespace", start)
        recipe = source[start:end]
        assert "term.fermion_order != std::vector<int>{0, 1, 2, 3}" in recipe, name
        assert recipe.count(
            "term_opts.orderExternalFermions = reorder_external_fermions;"
        ) >= 2, name
