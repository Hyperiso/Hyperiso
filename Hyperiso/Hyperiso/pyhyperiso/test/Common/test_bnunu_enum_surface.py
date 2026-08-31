"""Python enum-surface regression tests for B -> K(*) nu nu support."""

from pyhyperiso.Common import Decays, Observables, WCoeff, WGroup


def test_bnunu_public_enums_are_exposed():
    assert Decays.B__K_nu_nu.value is not None
    assert Decays.B__Kstar_nu_nu.value is not None

    assert Observables.BR_B__K_NU_NU.value is not None
    assert Observables.BR_B0__KS_NU_NU.value is not None
    assert Observables.BR_B__KSTAR_NU_NU.value is not None
    assert Observables.BR_B0__KSTAR0_NU_NU.value is not None

    assert WGroup.BNuNu.value is not None
    for name in (
        "CNU_L_EE", "CNU_R_EE", "CNU_L_EMU", "CNU_R_EMU",
        "CNU_L_ETAU", "CNU_R_ETAU", "CNU_L_MUE", "CNU_R_MUE",
        "CNU_L_MUMU", "CNU_R_MUMU", "CNU_L_MUTAU", "CNU_R_MUTAU",
        "CNU_L_TAUE", "CNU_R_TAUE", "CNU_L_TAUMU", "CNU_R_TAUMU",
        "CNU_L_TAUTAU", "CNU_R_TAUTAU",
    ):
        assert WCoeff[name].value is not None
