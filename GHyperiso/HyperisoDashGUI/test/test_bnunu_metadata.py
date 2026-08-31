"""GUI metadata checks for B -> K(*) nu nu support."""

from pyhyperiso_dash.latex_data.decay_name_to_latex_map import DECAY_ENUM_TO_LATEX_MAP
from pyhyperiso_dash.latex_data.observable_flha_to_latex_map import (
    OBSERVABLE_ENUM_TO_FLHA_MAP,
    OBSERVABLE_ENUM_TO_LATEX_MAP,
)
from pyhyperiso_dash.latex_data.wilson_lha_to_latex_map import (
    WILSON_ENUM_TO_LHA_ID_MAP,
    WILSON_NAME_TO_LATEX_MAP,
)


def test_bnunu_gui_labels_exist():
    assert "B__K_nu_nu" in DECAY_ENUM_TO_LATEX_MAP
    assert "B__Kstar_nu_nu" in DECAY_ENUM_TO_LATEX_MAP
    assert OBSERVABLE_ENUM_TO_FLHA_MAP["BR_B__K_NU_NU"] == (521, 1, 3, 321, 14, -14)
    assert OBSERVABLE_ENUM_TO_FLHA_MAP["BR_B0__KS_NU_NU"] == (511, 1, 3, 310, 14, -14)
    for name in (
        "BR_B__K_NU_NU",
        "BR_B0__KS_NU_NU",
        "BR_B__KSTAR_NU_NU",
        "BR_B0__KSTAR0_NU_NU",
    ):
        assert name in OBSERVABLE_ENUM_TO_LATEX_MAP
    assert "CNU_L_MUMU" in WILSON_NAME_TO_LATEX_MAP
    assert "CNU_R_MUMU" in WILSON_NAME_TO_LATEX_MAP
    assert WILSON_ENUM_TO_LHA_ID_MAP["CNU_L_MUMU"] == (3051414, 4141)
