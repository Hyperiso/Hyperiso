from __future__ import annotations

from importlib.resources import files


def _template_text(name: str) -> str:
    template = files("pyhyperiso").joinpath("assets", "template", "MARTY", name)
    return template.read_text(encoding="utf-8")


def test_packaged_c9_template_keeps_split_reg_prop_abi() -> None:
    source = _template_text("C9.cpp")

    assert "HYPERISO_MARTY_TEMPLATE_ABI: semileptonic-c9-tree-first-split-regprop" in source
    assert "HyperisoMartyC9LinkerSelection::NonPhotonVector" in source
    assert "HyperisoMartyC9LinkerSelection::PhotonOnly" in source
    assert "DiagramParticleType::External" in source
    assert "DiagramParticleType::Mediator" in source


def test_packaged_semileptonic_templates_keep_operator_normalization_abi() -> None:
    for name in ("C9.cpp", "C10.cpp", "CP9.cpp", "CP10.cpp"):
        source = _template_text(name)
        assert "HYPERISO_MARTY_OPERATOR_NORM_ABI" in source, name


def test_packaged_c10_template_is_tree_first() -> None:
    source = _template_text("C10.cpp")

    assert "semileptonic-c10-tree-first-full-4f-v9" in source
    assert "if (C10_tree == CSL_0)" in source
    assert "mty::Order::TreeLevel" in source
    assert "mty::Order::OneLoop" in source
