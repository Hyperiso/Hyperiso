import pytest
from pathlib import Path

from pyhyperiso.core.Common.GeneralEnum import Model
from pyhyperiso.core.Core.Config import ExternalFlag, PyConfig


def test_py_config_default_to_cpp():
    cfg = PyConfig()
    cpp_cfg = cfg.to_cpp()

    assert isinstance(cpp_cfg.flags, dict)
    for flag in ExternalFlag:
        assert cpp_cfg.flags[flag.value] is False

    assert cpp_cfg.model == Model.SM.value
    assert cpp_cfg.mty_model_name is None
    assert cpp_cfg.mty_model_path is None


def test_py_config_custom_values():
    cfg = PyConfig(
        flags={
            ExternalFlag.IS_LHA_SPECTRUM: True,
            ExternalFlag.HAS_WILSON_INPUT: False,
            ExternalFlag.HAS_TH_OBSERVABLE_INPUT: True,
            ExternalFlag.USE_MARTY: True,
        },
        model=Model.THDM,
        mty_model_name="THDM_Model",
        mty_model_path=Path("/tmp/THDM_Model.h")
    )

    cpp_cfg = cfg.to_cpp()

    assert cpp_cfg.flags[ExternalFlag.IS_LHA_SPECTRUM.value] is True
    assert cpp_cfg.flags[ExternalFlag.HAS_WILSON_INPUT.value] is False
    assert cpp_cfg.flags[ExternalFlag.HAS_TH_OBSERVABLE_INPUT.value] is True
    assert cpp_cfg.flags[ExternalFlag.USE_MARTY.value] is True

    assert cpp_cfg.model == Model.THDM.value
    assert cpp_cfg.mty_model_name == "THDM_Model"
    assert str(cpp_cfg.mty_model_path) == "/tmp/THDM_Model.h"
