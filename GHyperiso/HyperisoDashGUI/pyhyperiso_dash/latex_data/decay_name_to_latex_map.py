#!/usr/bin/env python3
"""Decay-name -> LaTeX map.

Generated from the Decays section of Map(4).cpp.

Main public mapping:
    DECAY_NAME_TO_LATEX_MAP[name] -> LaTeX label
"""

from __future__ import annotations

from typing import Dict, Optional


DECAY_NAME_TO_LATEX_MAP: Dict[str, str] = {
    'B__D_l_nu': '$B \\to D\\,\\ell\\,\\nu$',
    'B__Dstar_l_nu': '$B \\to D^{*}\\,\\ell\\,\\nu$',
    'B__K*_l_l': '$B \\to K^{*}\\,\\ell^{+}\\ell^{-}$',
    'B__K_l_l': '$B \\to K\\,\\ell^{+}\\ell^{-}$',
    'B__Kstar_gamma': '$B \\to K^{*}\\,\\gamma$',
    'B__Xs': '$B \\to X_s\\,\\gamma$',
    'B__Xs_ll': '$B \\to X_s\\,\\ell^{+}\\ell^{-}$',
    'B__l_l': '$B \\to \\ell^{+}\\ell^{-}$',
    'B__l_nu': '$B \\to \\ell\\,\\nu$',
    'Bs__phi_l_l': '$B_s \\to \\phi\\,\\ell^{+}\\ell^{-}$',
    'D__l_nu': '$D \\to \\ell\\,\\nu$',
    'Ds__l_nu': '$D_s \\to \\ell\\,\\nu$',
    'K__l_l': '$K \\to \\ell^{+}\\ell^{-}$',
    'K__l_nu': '$K \\to \\ell\\,\\nu$',
    'K__pi_nu_nu': '$K \\to \\pi\\,\\nu\\bar{\\nu}$',
    'Lambda_b__Lambda_l_l': '$\\Lambda_b \\to \\Lambda\\,\\ell^{+}\\ell^{-}$',
    'M0_Mix': '$M^0\\text{--}\\bar{M}^0\\;\\mathrm{mixing}$',
}
DECAY_ENUM_TO_NAME_MAP: Dict[str, str] = {
    'B__D_l_nu': 'B__D_l_nu',
    'B__Dstar_l_nu': 'B__Dstar_l_nu',
    'B__K_l_l': 'B__K_l_l',
    'B__Kstar_gamma': 'B__Kstar_gamma',
    'B__Kstar_l_l': 'B__K*_l_l',
    'B__Xs_gamma': 'B__Xs',
    'B__Xs_l_l': 'B__Xs_ll',
    'B__l_l': 'B__l_l',
    'B__l_nu': 'B__l_nu',
    'Bs__phi_l_l': 'Bs__phi_l_l',
    'D__l_nu': 'D__l_nu',
    'Ds__l_nu': 'Ds__l_nu',
    'K__l_l': 'K__l_l',
    'K__l_nu': 'K__l_nu',
    'K__pi_nu_nu': 'K__pi_nu_nu',
    'Lambda_b__Lambda_l_l': 'Lambda_b__Lambda_l_l',
    'M0_Mix': 'M0_Mix',
}
DECAY_ENUM_TO_LATEX_MAP: Dict[str, str] = {
    'B__D_l_nu': '$B \\to D\\,\\ell\\,\\nu$',
    'B__Dstar_l_nu': '$B \\to D^{*}\\,\\ell\\,\\nu$',
    'B__K_l_l': '$B \\to K\\,\\ell^{+}\\ell^{-}$',
    'B__Kstar_gamma': '$B \\to K^{*}\\,\\gamma$',
    'B__Kstar_l_l': '$B \\to K^{*}\\,\\ell^{+}\\ell^{-}$',
    'B__Xs_gamma': '$B \\to X_s\\,\\gamma$',
    'B__Xs_l_l': '$B \\to X_s\\,\\ell^{+}\\ell^{-}$',
    'B__l_l': '$B \\to \\ell^{+}\\ell^{-}$',
    'B__l_nu': '$B \\to \\ell\\,\\nu$',
    'Bs__phi_l_l': '$B_s \\to \\phi\\,\\ell^{+}\\ell^{-}$',
    'D__l_nu': '$D \\to \\ell\\,\\nu$',
    'Ds__l_nu': '$D_s \\to \\ell\\,\\nu$',
    'K__l_l': '$K \\to \\ell^{+}\\ell^{-}$',
    'K__l_nu': '$K \\to \\ell\\,\\nu$',
    'K__pi_nu_nu': '$K \\to \\pi\\,\\nu\\bar{\\nu}$',
    'Lambda_b__Lambda_l_l': '$\\Lambda_b \\to \\Lambda\\,\\ell^{+}\\ell^{-}$',
    'M0_Mix': '$M^0\\text{--}\\bar{M}^0\\;\\mathrm{mixing}$',
}
def get_decay_latex_name(name: str) -> Optional[str]:
    """Return the LaTeX label associated with an internal decay name."""
    return DECAY_NAME_TO_LATEX_MAP.get(name)
