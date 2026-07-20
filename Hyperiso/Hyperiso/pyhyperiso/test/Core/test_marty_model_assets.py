from __future__ import annotations

import json
import re
from importlib.resources import files


def _asset_text(*parts: str) -> str:
    return files("pyhyperiso").joinpath("assets", *parts).read_text(encoding="utf-8")


def test_packaged_thdm_mapping_covers_model_parameters() -> None:
    model = _asset_text("input_files", "marty_model", "thdm.h")
    mapping = json.loads(_asset_text("input_files", "marty_mapping", "thdm.json"))

    declared_lambdas = set(re.findall(r'constant_s\("(lambda_[1-7])"\)', model))
    assert declared_lambdas == {f"lambda_{index}" for index in range(1, 6)}
    assert declared_lambdas <= set(mapping)

    # Keep the runtime mapping synchronized with the source mapping, including
    # lambda_6/lambda_7 for compatible extended THDM headers.
    for index in range(1, 8):
        entry = mapping[f"lambda_{index}"]
        assert entry["block"] == "MINPAR"
        assert str(entry["pdgCode"]) == str(10 + index)
