#!/usr/bin/env python3
"""
Inverse nuisance-name map: new-format (block, id) targets to LaTeX labels.

Generated from nuisance_map(1).py by inverting:
- LEGACY_TO_NEW_NUISANCE_MAP
- LEGACY_TO_NEW_MULTI_MAP
- COMPRESSED_OR_SHARED_MAP

Values are LaTeX strings, ready to be rendered in math mode.
"""

from __future__ import annotations

from typing import Dict, Tuple

NewTarget = Tuple[str, str]

NEW_TO_LATEX_NUISANCE_MAP: Dict[NewTarget, str] = {
    ('SMINPUTS', '3'): r"$\alpha_s(M_Z)$",  # alphas_MZ
    ('SMINPUTS', '5'): r"$m_b$",  # mass_b
    ('MASS', '4'): r"$m_c$",  # mass_c
    ('MASS', '3'): r"$m_s$",  # mass_s
    ('SMINPUTS', '6'): r"$m_t^{\mathrm{pole}}$",  # mass_top_pole
    ('MASS', '25'): r"$m_{h^0}$",  # mass_h0
    ('VCKMIN', '1'): r"$\lambda_{\mathrm{CKM}}$",  # CKM_lambda
    ('VCKMIN', '2'): r"$A_{\mathrm{CKM}}$",  # CKM_A
    ('VCKMIN', '3'): r"$\bar{\rho}$",  # CKM_rhobar
    ('VCKMIN', '4'): r"$\bar{\eta}$",  # CKM_etabar
    ('B_Xs', '7'): r"$\mu_c(b\to s\gamma)$",  # mu_c_bsg
    ('B_Xs', '2'): r"$\mathrm{BR}(B\to X_c\ell\nu)_{\mathrm{exp}}$",  # BR_BXclnu_exp
    ('B_Xs', '3'): r"$\mu_G^2(b\to s\gamma)$",  # mu_G2_bsg
    ('B_Xs', '4'): r"$\rho_D^3(b\to s\gamma)$",  # rho_D3_bsg
    ('B_Xs', '5'): r"$\rho_{LS}^3(b\to s\gamma)$",  # rho_LS3_bsg
    ('B_Xs', '10'): r"$\mathrm{BR}(B\to X_s\gamma)_{\mathrm{rand}}$",  # bsgamma_rand
    ('B_Xsll', '3_15_0'): r"$\mathrm{BR}(B\to X_s \tau^+\tau^-)_{\mathrm{low}\ q^2}^{\mathrm{rand}}$",  # BRBXstautau_lowq2_rand
    ('B_Xsll', '3_15_1'): r"$\mathrm{BR}(B\to X_s \tau^+\tau^-)_{\mathrm{high}\ q^2}^{\mathrm{rand}}$",  # BRBXstautau_highq2_rand
    ('B_Xsll', '3_15_2'): r"$\mathrm{BR}(B\to X_s \tau^+\tau^-)_{\mathrm{full}}^{\mathrm{rand}}$",  # BRBXstautau_full_rand
    ('B_Xsll', '3_13_0'): r"$\mathrm{BR}(B\to X_s \mu^+\mu^-)_{\mathrm{low}\ q^2}^{\mathrm{rand}}$",  # BRBXsmumu_lowq2_rand
    ('B_Xsll', '3_13_1'): r"$\mathrm{BR}(B\to X_s \mu^+\mu^-)_{\mathrm{high}\ q^2}^{\mathrm{rand}}$",  # BRBXsmumu_highq2_rand
    ('B_Xsll', '3_13_2'): r"$\mathrm{BR}(B\to X_s \mu^+\mu^-)_{\mathrm{full}}^{\mathrm{rand}}$",  # BRBXsmumu_full_rand
    ('B_Xsll', '3_11_0'): r"$\mathrm{BR}(B\to X_s e^+e^-)_{\mathrm{low}\ q^2}^{\mathrm{rand}}$",  # BRBXsee_lowq2_rand
    ('B_Xsll', '3_11_1'): r"$\mathrm{BR}(B\to X_s e^+e^-)_{\mathrm{high}\ q^2}^{\mathrm{rand}}$",  # BRBXsee_highq2_rand
    ('B_Xsll', '3_11_2'): r"$\mathrm{BR}(B\to X_s e^+e^-)_{\mathrm{full}}^{\mathrm{rand}}$",  # BRBXsee_full_rand
    ('FCONST', '323_1'): r"$f_{K^*}^{\parallel}$",  # f_Kstar_par
    ('FCONST', '323_2'): r"$f_{K^*}^{\perp}$",  # f_Kstar_perp
    ('B_Ks', '7_1'): r"$a_1^{\perp,K^*}$",  # a1perp
    ('B_Ks', '7_2'): r"$a_2^{\perp,K^*}$",  # a2perp
    ('B_Ks', '8_1'): r"$a_1^{\parallel,K^*}$",  # a1par
    ('B_Ks', '8_2'): r"$a_2^{\parallel,K^*}$",  # a2par
    ('B_Ks', '16'): r"$T_1(B\to K^*)$",  # T1_BKstar
    ('B_Ks', '1_1_0'): r"$a_{0}^{A_0}(B\to K^*)$",  # a0A0_BKstar
    ('B_Ks', '1_1_1'): r"$a_{1}^{A_0}(B\to K^*)$",  # a1A0_BKstar
    ('B_Ks', '1_1_2'): r"$a_{2}^{A_0}(B\to K^*)$",  # a2A0_BKstar
    ('B_Ks', '1_2_0'): r"$a_{0}^{A_1}(B\to K^*)$",  # a0A1_BKstar
    ('B_Ks', '1_2_1'): r"$a_{1}^{A_1}(B\to K^*)$",  # a1A1_BKstar
    ('B_Ks', '1_2_2'): r"$a_{2}^{A_1}(B\to K^*)$",  # a2A1_BKstar
    ('B_Ks', '1_3_0'): r"$a_{0}^{A_{12}}(B\to K^*)$",  # a0A12_BKstar
    ('B_Ks', '1_3_1'): r"$a_{1}^{A_{12}}(B\to K^*)$",  # a1A12_BKstar
    ('B_Ks', '1_3_2'): r"$a_{2}^{A_{12}}(B\to K^*)$",  # a2A12_BKstar
    ('B_Ks', '1_4_0'): r"$a_{0}^{V}(B\to K^*)$",  # a0V_BKstar
    ('B_Ks', '1_4_1'): r"$a_{1}^{V}(B\to K^*)$",  # a1V_BKstar
    ('B_Ks', '1_4_2'): r"$a_{2}^{V}(B\to K^*)$",  # a2V_BKstar
    ('B_Ks', '1_5_0'): r"$a_{0}^{T_1}(B\to K^*)$",  # a0T1_BKstar
    ('B_Ks', '1_5_1'): r"$a_{1}^{T_1}(B\to K^*)$",  # a1T1_BKstar
    ('B_Ks', '1_5_2'): r"$a_{2}^{T_1}(B\to K^*)$",  # a2T1_BKstar
    ('B_Ks', '1_6_0'): r"$a_{0}^{T_2}(B\to K^*)$",  # a0T2_BKstar
    ('B_Ks', '1_6_1'): r"$a_{1}^{T_2}(B\to K^*)$",  # a1T2_BKstar
    ('B_Ks', '1_6_2'): r"$a_{2}^{T_2}(B\to K^*)$",  # a2T2_BKstar
    ('B_Ks', '1_7_0'): r"$a_{0}^{T_{23}}(B\to K^*)$",  # a0T23_BKstar
    ('B_Ks', '1_7_1'): r"$a_{1}^{T_{23}}(B\to K^*)$",  # a1T23_BKstar
    ('B_Ks', '1_7_2'): r"$a_{2}^{T_{23}}(B\to K^*)$",  # a2T23_BKstar
    ('B_Ks', '2_1_0'): r"$a_{0}^{A_0}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a0A0_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_1_1'): r"$a_{1}^{A_0}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a1A0_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_1_2'): r"$a_{2}^{A_0}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a2A0_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_2_0'): r"$a_{0}^{A_1}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a0A1_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_2_1'): r"$a_{1}^{A_1}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a1A1_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_2_2'): r"$a_{2}^{A_1}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a2A1_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_3_0'): r"$a_{0}^{A_{12}}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a0A12_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_3_1'): r"$a_{1}^{A_{12}}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a1A12_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_3_2'): r"$a_{2}^{A_{12}}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a2A12_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_4_0'): r"$a_{0}^{V}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a0V_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_4_1'): r"$a_{1}^{V}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a1V_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_4_2'): r"$a_{2}^{V}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a2V_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_5_0'): r"$a_{0}^{T_1}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a0T1_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_5_1'): r"$a_{1}^{T_1}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a1T1_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_5_2'): r"$a_{2}^{T_1}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a2T1_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_6_0'): r"$a_{0}^{T_2}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a0T2_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_6_1'): r"$a_{1}^{T_2}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a1T2_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_6_2'): r"$a_{2}^{T_2}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a2T2_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_7_0'): r"$a_{0}^{T_{23}}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a0T23_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_7_1'): r"$a_{1}^{T_{23}}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a1T23_BKstar_BSZ_LCSRonly
    ('B_Ks', '2_7_2'): r"$a_{2}^{T_{23}}(B\to K^*)_{\mathrm{BSZ,LCSR-only}}$",  # a2T23_BKstar_BSZ_LCSRonly
    ('B_Ks', '3_1_0'): r"$a_{0}^{A_0}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a0A0_BKstar_GRvDV_BSZ
    ('B_Ks', '3_1_1'): r"$a_{1}^{A_0}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a1A0_BKstar_GRvDV_BSZ
    ('B_Ks', '3_1_2'): r"$a_{2}^{A_0}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a2A0_BKstar_GRvDV_BSZ
    ('B_Ks', '3_2_0'): r"$a_{0}^{A_1}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a0A1_BKstar_GRvDV_BSZ
    ('B_Ks', '3_2_1'): r"$a_{1}^{A_1}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a1A1_BKstar_GRvDV_BSZ
    ('B_Ks', '3_2_2'): r"$a_{2}^{A_1}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a2A1_BKstar_GRvDV_BSZ
    ('B_Ks', '3_3_0'): r"$a_{0}^{A_{12}}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a0A12_BKstar_GRvDV_BSZ
    ('B_Ks', '3_3_1'): r"$a_{1}^{A_{12}}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a1A12_BKstar_GRvDV_BSZ
    ('B_Ks', '3_3_2'): r"$a_{2}^{A_{12}}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a2A12_BKstar_GRvDV_BSZ
    ('B_Ks', '3_4_0'): r"$a_{0}^{V}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a0V_BKstar_GRvDV_BSZ
    ('B_Ks', '3_4_1'): r"$a_{1}^{V}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a1V_BKstar_GRvDV_BSZ
    ('B_Ks', '3_4_2'): r"$a_{2}^{V}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a2V_BKstar_GRvDV_BSZ
    ('B_Ks', '3_5_0'): r"$a_{0}^{T_1}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a0T1_BKstar_GRvDV_BSZ
    ('B_Ks', '3_5_1'): r"$a_{1}^{T_1}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a1T1_BKstar_GRvDV_BSZ
    ('B_Ks', '3_5_2'): r"$a_{2}^{T_1}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a2T1_BKstar_GRvDV_BSZ
    ('B_Ks', '3_6_0'): r"$a_{0}^{T_2}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a0T2_BKstar_GRvDV_BSZ
    ('B_Ks', '3_6_1'): r"$a_{1}^{T_2}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a1T2_BKstar_GRvDV_BSZ
    ('B_Ks', '3_6_2'): r"$a_{2}^{T_2}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a2T2_BKstar_GRvDV_BSZ
    ('B_Ks', '3_7_0'): r"$a_{0}^{T_{23}}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a0T23_BKstar_GRvDV_BSZ
    ('B_Ks', '3_7_1'): r"$a_{1}^{T_{23}}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a1T23_BKstar_GRvDV_BSZ
    ('B_Ks', '3_7_2'): r"$a_{2}^{T_{23}}(B\to K^*)_{\mathrm{GRvDV,BSZ}}$",  # a2T23_BKstar_GRvDV_BSZ
    ('B_Ks', '4_1_0'): r"$a_{0}^{A_0}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0A0_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_1_1'): r"$a_{1}^{A_0}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1A0_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_1_2'): r"$a_{2}^{A_0}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2A0_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_2_0'): r"$a_{0}^{A_1}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0A1_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_2_1'): r"$a_{1}^{A_1}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1A1_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_2_2'): r"$a_{2}^{A_1}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2A1_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_3_0'): r"$a_{0}^{A_{12}}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0A12_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_3_1'): r"$a_{1}^{A_{12}}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1A12_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_3_2'): r"$a_{2}^{A_{12}}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2A12_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_4_0'): r"$a_{0}^{V}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0V_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_4_1'): r"$a_{1}^{V}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1V_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_4_2'): r"$a_{2}^{V}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2V_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_5_0'): r"$a_{0}^{T_1}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0T1_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_5_1'): r"$a_{1}^{T_1}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1T1_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_5_2'): r"$a_{2}^{T_1}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2T1_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_6_0'): r"$a_{0}^{T_2}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0T2_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_6_1'): r"$a_{1}^{T_2}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1T2_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_6_2'): r"$a_{2}^{T_2}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2T2_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_7_0'): r"$a_{0}^{T_{23}}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0T23_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_7_1'): r"$a_{1}^{T_{23}}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1T23_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '4_7_2'): r"$a_{2}^{T_{23}}(B\to K^*)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2T23_BKstar_GKvD_LCSR_Lattice
    ('B_Ks', '5_1_0'): r"$a_{0}^{A_0}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a0A0_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_1_1'): r"$a_{1}^{A_0}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a1A0_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_1_2'): r"$a_{2}^{A_0}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a2A0_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_2_0'): r"$a_{0}^{A_1}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a0A1_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_2_1'): r"$a_{1}^{A_1}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a1A1_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_2_2'): r"$a_{2}^{A_1}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a2A1_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_3_0'): r"$a_{0}^{A_{12}}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a0A12_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_3_1'): r"$a_{1}^{A_{12}}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a1A12_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_3_2'): r"$a_{2}^{A_{12}}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a2A12_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_4_0'): r"$a_{0}^{V}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a0V_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_4_1'): r"$a_{1}^{V}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a1V_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_4_2'): r"$a_{2}^{V}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a2V_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_5_0'): r"$a_{0}^{T_1}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a0T1_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_5_1'): r"$a_{1}^{T_1}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a1T1_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_5_2'): r"$a_{2}^{T_1}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a2T1_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_6_0'): r"$a_{0}^{T_2}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a0T2_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_6_1'): r"$a_{1}^{T_2}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a1T2_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_6_2'): r"$a_{2}^{T_2}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a2T2_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_7_0'): r"$a_{0}^{T_{23}}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a0T23_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_7_1'): r"$a_{1}^{T_{23}}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a1T23_BKstar_GKvD_LCSRonly
    ('B_Ks', '5_7_2'): r"$a_{2}^{T_{23}}(B\to K^*)_{\mathrm{GKvD,LCSR-only}}$",  # a2T23_BKstar_GKvD_LCSRonly
    ('B_Ks', '6'): r"$\sigma_{\mathrm{syst}}(B\to K^*)_{\mathrm{HLMW}}$",  # systErr_BKstar_HLMW
    ('B_Ks', '6_1_0'): r"$a_{0}^{A_0}(B\to K^*)_{\mathrm{HLMW}}$",  # a0A0_BKstar_HLMW
    ('B_Ks', '6_1_1'): r"$a_{1}^{A_0}(B\to K^*)_{\mathrm{HLMW}}$",  # a1A0_BKstar_HLMW
    ('B_Ks', '6_2_0'): r"$a_{0}^{A_1}(B\to K^*)_{\mathrm{HLMW}}$",  # a0A1_BKstar_HLMW
    ('B_Ks', '6_2_1'): r"$a_{1}^{A_1}(B\to K^*)_{\mathrm{HLMW}}$",  # a1A1_BKstar_HLMW
    ('B_Ks', '6_3_0'): r"$a_{0}^{A_{12}}(B\to K^*)_{\mathrm{HLMW}}$",  # a0A12_BKstar_HLMW
    ('B_Ks', '6_3_1'): r"$a_{1}^{A_{12}}(B\to K^*)_{\mathrm{HLMW}}$",  # a1A12_BKstar_HLMW
    ('B_Ks', '6_4_0'): r"$a_{0}^{V}(B\to K^*)_{\mathrm{HLMW}}$",  # a0V_BKstar_HLMW
    ('B_Ks', '6_4_1'): r"$a_{1}^{V}(B\to K^*)_{\mathrm{HLMW}}$",  # a1V_BKstar_HLMW
    ('B_Ks', '6_5_0'): r"$a_{0}^{T_1}(B\to K^*)_{\mathrm{HLMW}}$",  # a0T1_BKstar_HLMW
    ('B_Ks', '6_5_1'): r"$a_{1}^{T_1}(B\to K^*)_{\mathrm{HLMW}}$",  # a1T1_BKstar_HLMW
    ('B_Ks', '6_6_0'): r"$a_{0}^{T_2}(B\to K^*)_{\mathrm{HLMW}}$",  # a0T2_BKstar_HLMW
    ('B_Ks', '6_6_1'): r"$a_{1}^{T_2}(B\to K^*)_{\mathrm{HLMW}}$",  # a1T2_BKstar_HLMW
    ('B_Ks', '6_7_0'): r"$a_{0}^{T_{23}}(B\to K^*)_{\mathrm{HLMW}}$",  # a0T23_BKstar_HLMW
    ('B_Ks', '6_7_1'): r"$a_{1}^{T_{23}}(B\to K^*)_{\mathrm{HLMW}}$",  # a1T23_BKstar_HLMW
    ('B_Ks', '19_1_1_0'): r"$\mathrm{Re}(\alpha_{\perp,0})$",  # real_alpha_perp0
    ('B_Ks', '19_1_1_1'): r"$\mathrm{Re}(\alpha_{\perp,1})$",  # real_alpha_perp1
    ('B_Ks', '19_1_1_2'): r"$\mathrm{Re}(\alpha_{\perp,2})$",  # real_alpha_perp2
    ('B_Ks', '19_1_2_0'): r"$\mathrm{Re}(\alpha_{\parallel,0})$",  # real_alpha_par0
    ('B_Ks', '19_1_2_1'): r"$\mathrm{Re}(\alpha_{\parallel,1})$",  # real_alpha_par1
    ('B_Ks', '19_1_2_2'): r"$\mathrm{Re}(\alpha_{\parallel,2})$",  # real_alpha_par2
    ('B_Ks', '19_1_3_0'): r"$\mathrm{Re}(\alpha_{0,0})$",  # real_alpha_zero0
    ('B_Ks', '19_1_3_1'): r"$\mathrm{Re}(\alpha_{0,1})$",  # real_alpha_zero1
    ('B_Ks', '19_2_1_0'): r"$\mathrm{Im}(\alpha_{\perp,0})$",  # imag_alpha_perp0
    ('B_Ks', '19_2_1_1'): r"$\mathrm{Im}(\alpha_{\perp,1})$",  # imag_alpha_perp1
    ('B_Ks', '19_2_1_2'): r"$\mathrm{Im}(\alpha_{\perp,2})$",  # imag_alpha_perp2
    ('B_Ks', '19_2_2_0'): r"$\mathrm{Im}(\alpha_{\parallel,0})$",  # imag_alpha_par0
    ('B_Ks', '19_2_2_1'): r"$\mathrm{Im}(\alpha_{\parallel,1})$",  # imag_alpha_par1
    ('B_Ks', '19_2_2_2'): r"$\mathrm{Im}(\alpha_{\parallel,2})$",  # imag_alpha_par2
    ('B_Ks', '19_2_3_0'): r"$\mathrm{Im}(\alpha_{0,0})$",  # imag_alpha_zero0
    ('B_Ks', '19_2_3_1'): r"$\mathrm{Im}(\alpha_{0,1})$",  # imag_alpha_zero1
    ('B_Ks', '20_1'): r"$\Delta C_9^{M_1}(\bar q^2)$",  # DeltaC9_M1_q2bar
    ('B_Ks', '20_2'): r"$\Delta C_9^{M_2}(\bar q^2)$",  # DeltaC9_M2_q2bar
    ('B_Ks', '20_3'): r"$\Delta C_9^{M_3}(\bar q^2)$",  # DeltaC9_M3_q2bar
    ('B_Ks', '21_1_1'): r"$r_{1}^{M_1}$",  # r1_M1
    ('B_Ks', '21_1_2'): r"$r_{1}^{M_2}$",  # r1_M2
    ('B_Ks', '21_1_3'): r"$r_{1}^{M_3}$",  # r1_M3
    ('B_Ks', '21_2_1'): r"$r_{2}^{M_1}$",  # r2_M1
    ('B_Ks', '21_2_2'): r"$r_{2}^{M_2}$",  # r2_M2
    ('B_Ks', '21_2_3'): r"$r_{2}^{M_3}$",  # r2_M3
    ('B_K', '1_1_0'): r"$a_{0}^{f_+}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a0fp_BK_AS_LCSR_Lattice
    ('B_K', '1_1_1'): r"$a_{1}^{f_+}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a1fp_BK_AS_LCSR_Lattice
    ('B_K', '1_1_2'): r"$a_{2}^{f_+}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a2fp_BK_AS_LCSR_Lattice
    ('B_K', '1_2_0'): r"$a_{0}^{f_0}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a0f0_BK_AS_LCSR_Lattice
    ('B_K', '1_2_1'): r"$a_{1}^{f_0}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a1f0_BK_AS_LCSR_Lattice
    ('B_K', '1_2_2'): r"$a_{2}^{f_0}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a2f0_BK_AS_LCSR_Lattice
    ('B_K', '1_2_3'): r"$a_{3}^{f_0}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a3f0_BK_AS_LCSR_Lattice
    ('B_K', '1_3_0'): r"$a_{0}^{f_T}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a0fT_BK_AS_LCSR_Lattice
    ('B_K', '1_3_1'): r"$a_{1}^{f_T}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a1fT_BK_AS_LCSR_Lattice
    ('B_K', '1_3_2'): r"$a_{2}^{f_T}(B\to K)_{\mathrm{AS,LCSR+Lattice}}$",  # a2fT_BK_AS_LCSR_Lattice
    ('B_K', '2_1_0'): r"$a_{0}^{f_+}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a0fp_BK_GRvDV_BSZ
    ('B_K', '2_1_1'): r"$a_{1}^{f_+}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a1fp_BK_GRvDV_BSZ
    ('B_K', '2_1_2'): r"$a_{2}^{f_+}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a2fp_BK_GRvDV_BSZ
    ('B_K', '2_2_1'): r"$a_{1}^{f_0}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a1f0_BK_GRvDV_BSZ
    ('B_K', '2_2_2'): r"$a_{2}^{f_0}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a2f0_BK_GRvDV_BSZ
    ('B_K', '2_3_0'): r"$a_{0}^{f_T}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a0fT_BK_GRvDV_BSZ
    ('B_K', '2_3_1'): r"$a_{1}^{f_T}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a1fT_BK_GRvDV_BSZ
    ('B_K', '2_3_2'): r"$a_{2}^{f_T}(B\to K)_{\mathrm{GRvDV,BSZ}}$",  # a2fT_BK_GRvDV_BSZ
    ('B_K', '3_1_0'): r"$a_{0}^{f_+}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0fp_BK_GKvD_LCSR_Lattice
    ('B_K', '3_1_1'): r"$a_{1}^{f_+}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1fp_BK_GKvD_LCSR_Lattice
    ('B_K', '3_1_2'): r"$a_{2}^{f_+}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2fp_BK_GKvD_LCSR_Lattice
    ('B_K', '3_2_1'): r"$a_{1}^{f_0}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1f0_BK_GKvD_LCSR_Lattice
    ('B_K', '3_2_2'): r"$a_{2}^{f_0}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2f0_BK_GKvD_LCSR_Lattice
    ('B_K', '3_3_0'): r"$a_{0}^{f_T}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a0fT_BK_GKvD_LCSR_Lattice
    ('B_K', '3_3_1'): r"$a_{1}^{f_T}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a1fT_BK_GKvD_LCSR_Lattice
    ('B_K', '3_3_2'): r"$a_{2}^{f_T}(B\to K)_{\mathrm{GKvD,LCSR+Lattice}}$",  # a2fT_BK_GKvD_LCSR_Lattice
    ('B_K', '4_1_0'): r"$a_{0}^{f_+}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a0fp_BK_GKvD_LCSRonly
    ('B_K', '4_1_1'): r"$a_{1}^{f_+}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a1fp_BK_GKvD_LCSRonly
    ('B_K', '4_1_2'): r"$a_{2}^{f_+}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a2fp_BK_GKvD_LCSRonly
    ('B_K', '4_2_1'): r"$a_{1}^{f_0}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a1f0_BK_GKvD_LCSRonly
    ('B_K', '4_2_2'): r"$a_{2}^{f_0}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a2f0_BK_GKvD_LCSRonly
    ('B_K', '4_3_0'): r"$a_{0}^{f_T}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a0fT_BK_GKvD_LCSRonly
    ('B_K', '4_3_1'): r"$a_{1}^{f_T}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a1fT_BK_GKvD_LCSRonly
    ('B_K', '4_3_2'): r"$a_{2}^{f_T}(B\to K)_{\mathrm{GKvD,LCSR-only}}$",  # a2fT_BK_GKvD_LCSRonly
    ('B_K', '5_1_0'): r"$a_{0}^{f_+}(B\to K)_{\mathrm{FLAG24}}$",  # a0fp_BK_FLAG24
    ('B_K', '5_1_1'): r"$a_{1}^{f_+}(B\to K)_{\mathrm{FLAG24}}$",  # a1fp_BK_FLAG24
    ('B_K', '5_1_2'): r"$a_{2}^{f_+}(B\to K)_{\mathrm{FLAG24}}$",  # a2fp_BK_FLAG24
    ('B_K', '5_2_0'): r"$a_{0}^{f_0}(B\to K)_{\mathrm{FLAG24}}$",  # a0f0_BK_FLAG24
    ('B_K', '5_2_1'): r"$a_{1}^{f_0}(B\to K)_{\mathrm{FLAG24}}$",  # a1f0_BK_FLAG24
    ('B_K', '5_3_0'): r"$a_{0}^{f_T}(B\to K)_{\mathrm{FLAG24}}$",  # a0fT_BK_FLAG24
    ('B_K', '5_3_1'): r"$a_{1}^{f_T}(B\to K)_{\mathrm{FLAG24}}$",  # a1fT_BK_FLAG24
    ('B_K', '5_3_2'): r"$a_{2}^{f_T}(B\to K)_{\mathrm{FLAG24}}$",  # a2fT_BK_FLAG24
    ('B_K', '6_1_1'): r"$a_{1}^{f_+}(B\to K)_{\mathrm{HPQCD22}}$",  # a1fp_BK_HPQCD22
    ('B_K', '6_1_2'): r"$a_{2}^{f_+}(B\to K)_{\mathrm{HPQCD22}}$",  # a2fp_BK_HPQCD22
    ('B_K', '6_2_0'): r"$a_{0}^{f_0}(B\to K)_{\mathrm{HPQCD22}}$",  # a0f0_BK_HPQCD22
    ('B_K', '6_2_1'): r"$a_{1}^{f_0}(B\to K)_{\mathrm{HPQCD22}}$",  # a1f0_BK_HPQCD22
    ('B_K', '6_2_2'): r"$a_{2}^{f_0}(B\to K)_{\mathrm{HPQCD22}}$",  # a2f0_BK_HPQCD22
    ('B_K', '6_3_0'): r"$a_{0}^{f_T}(B\to K)_{\mathrm{HPQCD22}}$",  # a0fT_BK_HPQCD22
    ('B_K', '6_3_1'): r"$a_{1}^{f_T}(B\to K)_{\mathrm{HPQCD22}}$",  # a1fT_BK_HPQCD22
    ('B_K', '6_3_2'): r"$a_{2}^{f_T}(B\to K)_{\mathrm{HPQCD22}}$",  # a2fT_BK_HPQCD22
    ('B_K', '8_1'): r"$a_1^K$",  # a1K
    ('B_K', '8_2'): r"$a_2^K$",  # a2K
    ('FLIFE', '531'): r"$\tau_{B_s}$",  # life_Bs
    ('FCONST', '531_1'): r"$f_{B_s}$",  # f_Bs
    ('B_ll', '1'): r"$y_s$",  # ys_Bs
    ('FCONST', '333_1'): r"$f_{\phi}^{\parallel}$",  # f_phi_par
    ('FCONST', '333_2'): r"$f_{\phi}^{\perp}$",  # f_phi_perp
    ('B_phi', '7_1'): r"$a_1^{\perp,\phi}$",  # a1phi_perp
    ('B_phi', '7_2'): r"$a_2^{\perp,\phi}$",  # a2phi_perp
    ('B_phi', '8_1'): r"$a_1^{\parallel,\phi}$",  # a1phi_par
    ('B_phi', '8_2'): r"$a_2^{\parallel,\phi}$",  # a2phi_par
    ('B_phi', '1_1_0'): r"$a_{0}^{A_0}(B_s\to\phi)$",  # a0A0_Bsphi
    ('B_phi', '1_1_1'): r"$a_{1}^{A_0}(B_s\to\phi)$",  # a1A0_Bsphi
    ('B_phi', '1_1_2'): r"$a_{2}^{A_0}(B_s\to\phi)$",  # a2A0_Bsphi
    ('B_phi', '1_2_0'): r"$a_{0}^{A_1}(B_s\to\phi)$",  # a0A1_Bsphi
    ('B_phi', '1_2_1'): r"$a_{1}^{A_1}(B_s\to\phi)$",  # a1A1_Bsphi
    ('B_phi', '1_2_2'): r"$a_{2}^{A_1}(B_s\to\phi)$",  # a2A1_Bsphi
    ('B_phi', '1_3_0'): r"$a_{0}^{A_{12}}(B_s\to\phi)$",  # a0A12_Bsphi
    ('B_phi', '1_3_1'): r"$a_{1}^{A_{12}}(B_s\to\phi)$",  # a1A12_Bsphi
    ('B_phi', '1_3_2'): r"$a_{2}^{A_{12}}(B_s\to\phi)$",  # a2A12_Bsphi
    ('B_phi', '1_4_0'): r"$a_{0}^{V}(B_s\to\phi)$",  # a0V_Bsphi
    ('B_phi', '1_4_1'): r"$a_{1}^{V}(B_s\to\phi)$",  # a1V_Bsphi
    ('B_phi', '1_4_2'): r"$a_{2}^{V}(B_s\to\phi)$",  # a2V_Bsphi
    ('B_phi', '1_5_0'): r"$a_{0}^{T_1}(B_s\to\phi)$",  # a0T1_Bsphi
    ('B_phi', '1_5_1'): r"$a_{1}^{T_1}(B_s\to\phi)$",  # a1T1_Bsphi
    ('B_phi', '1_5_2'): r"$a_{2}^{T_1}(B_s\to\phi)$",  # a2T1_Bsphi
    ('B_phi', '1_6_0'): r"$a_{0}^{T_2}(B_s\to\phi)$",  # a0T2_Bsphi
    ('B_phi', '1_6_1'): r"$a_{1}^{T_2}(B_s\to\phi)$",  # a1T2_Bsphi
    ('B_phi', '1_6_2'): r"$a_{2}^{T_2}(B_s\to\phi)$",  # a2T2_Bsphi
    ('B_phi', '1_7_0'): r"$a_{0}^{T_{23}}(B_s\to\phi)$",  # a0T23_Bsphi
    ('B_phi', '1_7_1'): r"$a_{1}^{T_{23}}(B_s\to\phi)$",  # a1T23_Bsphi
    ('B_phi', '1_7_2'): r"$a_{2}^{T_{23}}(B_s\to\phi)$",  # a2T23_Bsphi
    ('B_phi', '2_1_0'): r"$a_{0}^{A_0}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a0A0_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_1_1'): r"$a_{1}^{A_0}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a1A0_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_1_2'): r"$a_{2}^{A_0}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a2A0_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_2_0'): r"$a_{0}^{A_1}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a0A1_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_2_1'): r"$a_{1}^{A_1}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a1A1_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_2_2'): r"$a_{2}^{A_1}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a2A1_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_3_0'): r"$a_{0}^{A_{12}}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a0A12_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_3_1'): r"$a_{1}^{A_{12}}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a1A12_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_3_2'): r"$a_{2}^{A_{12}}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a2A12_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_4_0'): r"$a_{0}^{V}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a0V_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_4_1'): r"$a_{1}^{V}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a1V_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_4_2'): r"$a_{2}^{V}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a2V_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_5_0'): r"$a_{0}^{T_1}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a0T1_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_5_1'): r"$a_{1}^{T_1}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a1T1_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_5_2'): r"$a_{2}^{T_1}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a2T1_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_6_0'): r"$a_{0}^{T_2}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a0T2_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_6_1'): r"$a_{1}^{T_2}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a1T2_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_6_2'): r"$a_{2}^{T_2}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a2T2_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_7_0'): r"$a_{0}^{T_{23}}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a0T23_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_7_1'): r"$a_{1}^{T_{23}}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a1T23_Bsphi_BSZ_LCSRonly
    ('B_phi', '2_7_2'): r"$a_{2}^{T_{23}}(B_s\to\phi)_{\mathrm{BSZ,LCSR-only}}$",  # a2T23_Bsphi_BSZ_LCSRonly
    ('B_phi', '3_1_0'): r"$a_{0}^{A_0}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a0A0_Bsphi_GRvDV_BSZ
    ('B_phi', '3_1_1'): r"$a_{1}^{A_0}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a1A0_Bsphi_GRvDV_BSZ
    ('B_phi', '3_1_2'): r"$a_{2}^{A_0}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a2A0_Bsphi_GRvDV_BSZ
    ('B_phi', '3_2_0'): r"$a_{0}^{A_1}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a0A1_Bsphi_GRvDV_BSZ
    ('B_phi', '3_2_1'): r"$a_{1}^{A_1}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a1A1_Bsphi_GRvDV_BSZ
    ('B_phi', '3_2_2'): r"$a_{2}^{A_1}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a2A1_Bsphi_GRvDV_BSZ
    ('B_phi', '3_3_0'): r"$a_{0}^{A_{12}}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a0A12_Bsphi_GRvDV_BSZ
    ('B_phi', '3_3_1'): r"$a_{1}^{A_{12}}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a1A12_Bsphi_GRvDV_BSZ
    ('B_phi', '3_3_2'): r"$a_{2}^{A_{12}}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a2A12_Bsphi_GRvDV_BSZ
    ('B_phi', '3_4_0'): r"$a_{0}^{V}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a0V_Bsphi_GRvDV_BSZ
    ('B_phi', '3_4_1'): r"$a_{1}^{V}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a1V_Bsphi_GRvDV_BSZ
    ('B_phi', '3_4_2'): r"$a_{2}^{V}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a2V_Bsphi_GRvDV_BSZ
    ('B_phi', '3_5_0'): r"$a_{0}^{T_1}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a0T1_Bsphi_GRvDV_BSZ
    ('B_phi', '3_5_1'): r"$a_{1}^{T_1}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a1T1_Bsphi_GRvDV_BSZ
    ('B_phi', '3_5_2'): r"$a_{2}^{T_1}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a2T1_Bsphi_GRvDV_BSZ
    ('B_phi', '3_6_0'): r"$a_{0}^{T_2}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a0T2_Bsphi_GRvDV_BSZ
    ('B_phi', '3_6_1'): r"$a_{1}^{T_2}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a1T2_Bsphi_GRvDV_BSZ
    ('B_phi', '3_6_2'): r"$a_{2}^{T_2}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a2T2_Bsphi_GRvDV_BSZ
    ('B_phi', '3_7_0'): r"$a_{0}^{T_{23}}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a0T23_Bsphi_GRvDV_BSZ
    ('B_phi', '3_7_1'): r"$a_{1}^{T_{23}}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a1T23_Bsphi_GRvDV_BSZ
    ('B_phi', '3_7_2'): r"$a_{2}^{T_{23}}(B_s\to\phi)_{\mathrm{GRvDV,BSZ}}$",  # a2T23_Bsphi_GRvDV_BSZ
    ('B_phi', '6'): r"$\sigma_{\mathrm{syst}}(B_s\to\phi)_{\mathrm{HLMW}}$",  # systErr_Bsphi_HLMW
    ('B_phi', '6_1_0'): r"$a_{0}^{A_0}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a0A0_Bsphi_HLMW
    ('B_phi', '6_1_1'): r"$a_{1}^{A_0}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a1A0_Bsphi_HLMW
    ('B_phi', '6_2_0'): r"$a_{0}^{A_1}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a0A1_Bsphi_HLMW
    ('B_phi', '6_2_1'): r"$a_{1}^{A_1}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a1A1_Bsphi_HLMW
    ('B_phi', '6_3_0'): r"$a_{0}^{A_{12}}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a0A12_Bsphi_HLMW
    ('B_phi', '6_3_1'): r"$a_{1}^{A_{12}}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a1A12_Bsphi_HLMW
    ('B_phi', '6_4_0'): r"$a_{0}^{V}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a0V_Bsphi_HLMW
    ('B_phi', '6_4_1'): r"$a_{1}^{V}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a1V_Bsphi_HLMW
    ('B_phi', '6_5_0'): r"$a_{0}^{T_1}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a0T1_Bsphi_HLMW
    ('B_phi', '6_5_1'): r"$a_{1}^{T_1}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a1T1_Bsphi_HLMW
    ('B_phi', '6_6_0'): r"$a_{0}^{T_2}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a0T2_Bsphi_HLMW
    ('B_phi', '6_6_1'): r"$a_{1}^{T_2}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a1T2_Bsphi_HLMW
    ('B_phi', '6_7_0'): r"$a_{0}^{T_{23}}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a0T23_Bsphi_HLMW
    ('B_phi', '6_7_1'): r"$a_{1}^{T_{23}}(B_s\to\phi)_{\mathrm{HLMW}}$",  # a1T23_Bsphi_HLMW
    ('FLIFE', '5122'): r"$\tau_{\Lambda_b}$",  # life_Lb
    ('Lb_L', '8'): r"$\alpha_{\Lambda}(\Lambda_b\to\Lambda\ell\ell)$",  # alphaL_LbLll
    ('Lb_L', '1_1_0'): r"$a_{0}^{f_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_fperp_LbLll
    ('Lb_L', '1_1_1'): r"$a_{1}^{f_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_fperp_LbLll
    ('Lb_L', '1_1_2'): r"$a_{2}^{f_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_fperp_LbLll
    ('Lb_L', '1_2_0'): r"$a_{0}^{h_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_hperp_LbLll
    ('Lb_L', '1_2_1'): r"$a_{1}^{h_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_hperp_LbLll
    ('Lb_L', '1_2_2'): r"$a_{2}^{h_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_hperp_LbLll
    ('Lb_L', '1_3_1'): r"$a_{1}^{g_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_gperp_LbLll
    ('Lb_L', '1_3_2'): r"$a_{2}^{g_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_gperp_LbLll
    ('Lb_L', '1_4_0'): r"$a_{0}^{f_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_fplus_LbLll
    ('Lb_L', '1_4_1'): r"$a_{1}^{f_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_fplus_LbLll
    ('Lb_L', '1_4_2'): r"$a_{2}^{f_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_fplus_LbLll
    ('Lb_L', '1_5_0'): r"$a_{0}^{h_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_hplus_LbLll
    ('Lb_L', '1_5_1'): r"$a_{1}^{h_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_hplus_LbLll
    ('Lb_L', '1_5_2'): r"$a_{2}^{h_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_hplus_LbLll
    ('Lb_L', '1_6_1'): r"$a_{1}^{g_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_gplus_LbLll
    ('Lb_L', '1_6_2'): r"$a_{2}^{g_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_gplus_LbLll
    ('Lb_L', '1_7_1'): r"$a_{1}^{\tilde h_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_htildeperp_LbLll
    ('Lb_L', '1_7_2'): r"$a_{2}^{\tilde h_\perp}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_htildeperp_LbLll
    ('Lb_L', '1_8_1'): r"$a_{1}^{\tilde h_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a1_HO_htildeplus_LbLll
    ('Lb_L', '1_8_2'): r"$a_{2}^{\tilde h_+}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a2_HO_htildeplus_LbLll
    ('K_ll', '1'): r"$\mathrm{BR}(K_L\to\gamma\gamma)_{\mathrm{exp}}$",  # BR_KLgammagamma_exp
    ('K_ll', '2'): r"$\mathrm{BR}(K_S\to\gamma\gamma)_{\mathrm{exp}}$",  # BR_KSgammagamma_exp
    ('B_Ks', '18_1_1'): r"$\delta A_L^\perp(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_ALperp_err_noq2
    ('B_Ks', '18_1_2'): r"$\delta A_R^\perp(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_ARperp_err_noq2
    ('B_Ks', '18_1_3'): r"$\delta A_L^\parallel(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_ALpar_err_noq2
    ('B_Ks', '18_1_4'): r"$\delta A_R^\parallel(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_ARpar_err_noq2
    ('B_Ks', '18_1_5'): r"$\delta A_L^0(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_AL0_err_noq2
    ('B_Ks', '18_1_6'): r"$\delta A_R^0(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_AR0_err_noq2
    ('B_Ks', '18_1_7'): r"$\delta A_t(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_At_err_noq2
    ('B_Ks', '18_1_8'): r"$\delta A_S(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKstarlow_AS_err_noq2
    ('B_Ks', '18_2_1'): r"$\delta A_L^\perp(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_ALperp_err_q2
    ('B_Ks', '18_2_2'): r"$\delta A_R^\perp(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_ARperp_err_q2
    ('B_Ks', '18_2_3'): r"$\delta A_L^\parallel(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_ALpar_err_q2
    ('B_Ks', '18_2_4'): r"$\delta A_R^\parallel(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_ARpar_err_q2
    ('B_Ks', '18_2_5'): r"$\delta A_L^0(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_AL0_err_q2
    ('B_Ks', '18_2_6'): r"$\delta A_R^0(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_AR0_err_q2
    ('B_Ks', '18_2_7'): r"$\delta A_t(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_At_err_q2
    ('B_Ks', '18_2_8'): r"$\delta A_S(B\to K^*)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKstarlow_AS_err_q2
    ('B_Ks', '18_3_1'): r"$\delta A_L^\perp(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_ALperp_err
    ('B_Ks', '18_3_2'): r"$\delta A_R^\perp(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_ARperp_err
    ('B_Ks', '18_3_3'): r"$\delta A_L^\parallel(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_ALpar_err
    ('B_Ks', '18_3_4'): r"$\delta A_R^\parallel(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_ARpar_err
    ('B_Ks', '18_3_5'): r"$\delta A_L^0(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_AL0_err
    ('B_Ks', '18_3_6'): r"$\delta A_R^0(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_AR0_err
    ('B_Ks', '18_3_7'): r"$\delta A_t(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_At_err
    ('B_Ks', '18_3_8'): r"$\delta A_S(B\to K^*)_{\mathrm{high}}$",  # BtoKstarhigh_AS_err
    ('B_K', '18_1_1'): r"$\delta F_V(B\to K)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKlow_FV_err_noq2
    ('B_K', '18_1_2'): r"$\delta F_A(B\to K)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKlow_FA_err_noq2
    ('B_K', '18_1_3'): r"$\delta F_S(B\to K)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKlow_FS_err_noq2
    ('B_K', '18_1_4'): r"$\delta F_P(B\to K)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # BtoKlow_FP_err_noq2
    ('B_K', '18_2_1'): r"$\delta F_V(B\to K)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKlow_FV_err_q2
    ('B_K', '18_2_2'): r"$\delta F_A(B\to K)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKlow_FA_err_q2
    ('B_K', '18_2_3'): r"$\delta F_S(B\to K)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKlow_FS_err_q2
    ('B_K', '18_2_4'): r"$\delta F_P(B\to K)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # BtoKlow_FP_err_q2
    ('B_K', '18_3_1'): r"$\delta F_V(B\to K)_{\mathrm{high}}$",  # BtoKhigh_FV_err
    ('B_K', '18_3_2'): r"$\delta F_A(B\to K)_{\mathrm{high}}$",  # BtoKhigh_FA_err
    ('B_K', '18_3_3'): r"$\delta F_S(B\to K)_{\mathrm{high}}$",  # BtoKhigh_FS_err
    ('B_K', '18_3_4'): r"$\delta F_P(B\to K)_{\mathrm{high}}$",  # BtoKhigh_FP_err
    ('B_phi', '18_1_1'): r"$\delta A_L^\perp(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_ALperp_err_noq2
    ('B_phi', '18_1_2'): r"$\delta A_R^\perp(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_ARperp_err_noq2
    ('B_phi', '18_1_3'): r"$\delta A_L^\parallel(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_ALpar_err_noq2
    ('B_phi', '18_1_4'): r"$\delta A_R^\parallel(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_ARpar_err_noq2
    ('B_phi', '18_1_5'): r"$\delta A_L^0(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_AL0_err_noq2
    ('B_phi', '18_1_6'): r"$\delta A_R^0(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_AR0_err_noq2
    ('B_phi', '18_1_7'): r"$\delta A_t(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_At_err_noq2
    ('B_phi', '18_1_8'): r"$\delta A_S(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-indep}}$",  # Bstophilow_AS_err_noq2
    ('B_phi', '18_2_1'): r"$\delta A_L^\perp(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_ALperp_err_q2
    ('B_phi', '18_2_2'): r"$\delta A_R^\perp(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_ARperp_err_q2
    ('B_phi', '18_2_3'): r"$\delta A_L^\parallel(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_ALpar_err_q2
    ('B_phi', '18_2_4'): r"$\delta A_R^\parallel(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_ARpar_err_q2
    ('B_phi', '18_2_5'): r"$\delta A_L^0(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_AL0_err_q2
    ('B_phi', '18_2_6'): r"$\delta A_R^0(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_AR0_err_q2
    ('B_phi', '18_2_7'): r"$\delta A_t(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_At_err_q2
    ('B_phi', '18_2_8'): r"$\delta A_S(B_s\to\phi)_{\mathrm{low},\ q^2\mathrm{-dep}}$",  # Bstophilow_AS_err_q2
    ('B_phi', '18_3_1'): r"$\delta A_L^\perp(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_ALperp_err
    ('B_phi', '18_3_2'): r"$\delta A_R^\perp(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_ARperp_err
    ('B_phi', '18_3_3'): r"$\delta A_L^\parallel(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_ALpar_err
    ('B_phi', '18_3_4'): r"$\delta A_R^\parallel(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_ARpar_err
    ('B_phi', '18_3_5'): r"$\delta A_L^0(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_AL0_err
    ('B_phi', '18_3_6'): r"$\delta A_R^0(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_AR0_err
    ('B_phi', '18_3_7'): r"$\delta A_t(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_At_err
    ('B_phi', '18_3_8'): r"$\delta A_S(B_s\to\phi)_{\mathrm{high}}$",  # Bstophihigh_AS_err
    ('FCONST', '511_1'): r"$f_B$",  # f_B
    ('FCONST', '521_1'): r"$f_B$",  # f_B
    ('FCONST', '311_1'): r"$f_K$",  # f_K
    ('FCONST', '321_1'): r"$f_K$",  # f_K
    ('B_K', '13'): r"$\lambda_{B^+}$",  # lambda_Bp
    ('B_Ks', '13'): r"$\lambda_{B^+}$",  # lambda_Bp
    ('B_phi', '13'): r"$\lambda_{B_s}$",  # lambda_Bsp
    ('Lb_L', '1_3_0'): r"$a_{0}^{g_{\mathrm{pp}}}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_gpp_LbLll
    ('Lb_L', '1_6_0'): r"$a_{0}^{g_{\mathrm{pp}}}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_gpp_LbLll
    ('Lb_L', '1_7_0'): r"$a_{0}^{\tilde h_{\mathrm{pp}}}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_htildepp_LbLll
    ('Lb_L', '1_8_0'): r"$a_{0}^{\tilde h_{\mathrm{pp}}}(\Lambda_b\to\Lambda\ell\ell)_{\mathrm{HO}}$",  # a0_HO_htildepp_LbLll
}

NEW_TO_LEGACY_NUISANCE_MAP: Dict[NewTarget, str] = {
    ('SMINPUTS', '3'): 'alphas_MZ',
    ('SMINPUTS', '5'): 'mass_b',
    ('MASS', '4'): 'mass_c',
    ('MASS', '3'): 'mass_s',
    ('SMINPUTS', '6'): 'mass_top_pole',
    ('MASS', '25'): 'mass_h0',
    ('VCKMIN', '1'): 'CKM_lambda',
    ('VCKMIN', '2'): 'CKM_A',
    ('VCKMIN', '3'): 'CKM_rhobar',
    ('VCKMIN', '4'): 'CKM_etabar',
    ('B_Xs', '7'): 'mu_c_bsg',
    ('B_Xs', '2'): 'BR_BXclnu_exp',
    ('B_Xs', '3'): 'mu_G2_bsg',
    ('B_Xs', '4'): 'rho_D3_bsg',
    ('B_Xs', '5'): 'rho_LS3_bsg',
    ('B_Xs', '10'): 'bsgamma_rand',
    ('B_Xsll', '3_15_0'): 'BRBXstautau_lowq2_rand',
    ('B_Xsll', '3_15_1'): 'BRBXstautau_highq2_rand',
    ('B_Xsll', '3_15_2'): 'BRBXstautau_full_rand',
    ('B_Xsll', '3_13_0'): 'BRBXsmumu_lowq2_rand',
    ('B_Xsll', '3_13_1'): 'BRBXsmumu_highq2_rand',
    ('B_Xsll', '3_13_2'): 'BRBXsmumu_full_rand',
    ('B_Xsll', '3_11_0'): 'BRBXsee_lowq2_rand',
    ('B_Xsll', '3_11_1'): 'BRBXsee_highq2_rand',
    ('B_Xsll', '3_11_2'): 'BRBXsee_full_rand',
    ('FCONST', '323_1'): 'f_Kstar_par',
    ('FCONST', '323_2'): 'f_Kstar_perp',
    ('B_Ks', '7_1'): 'a1perp',
    ('B_Ks', '7_2'): 'a2perp',
    ('B_Ks', '8_1'): 'a1par',
    ('B_Ks', '8_2'): 'a2par',
    ('B_Ks', '16'): 'T1_BKstar',
    ('B_Ks', '1_1_0'): 'a0A0_BKstar',
    ('B_Ks', '1_1_1'): 'a1A0_BKstar',
    ('B_Ks', '1_1_2'): 'a2A0_BKstar',
    ('B_Ks', '1_2_0'): 'a0A1_BKstar',
    ('B_Ks', '1_2_1'): 'a1A1_BKstar',
    ('B_Ks', '1_2_2'): 'a2A1_BKstar',
    ('B_Ks', '1_3_0'): 'a0A12_BKstar',
    ('B_Ks', '1_3_1'): 'a1A12_BKstar',
    ('B_Ks', '1_3_2'): 'a2A12_BKstar',
    ('B_Ks', '1_4_0'): 'a0V_BKstar',
    ('B_Ks', '1_4_1'): 'a1V_BKstar',
    ('B_Ks', '1_4_2'): 'a2V_BKstar',
    ('B_Ks', '1_5_0'): 'a0T1_BKstar',
    ('B_Ks', '1_5_1'): 'a1T1_BKstar',
    ('B_Ks', '1_5_2'): 'a2T1_BKstar',
    ('B_Ks', '1_6_0'): 'a0T2_BKstar',
    ('B_Ks', '1_6_1'): 'a1T2_BKstar',
    ('B_Ks', '1_6_2'): 'a2T2_BKstar',
    ('B_Ks', '1_7_0'): 'a0T23_BKstar',
    ('B_Ks', '1_7_1'): 'a1T23_BKstar',
    ('B_Ks', '1_7_2'): 'a2T23_BKstar',
    ('B_Ks', '2_1_0'): 'a0A0_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_1_1'): 'a1A0_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_1_2'): 'a2A0_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_2_0'): 'a0A1_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_2_1'): 'a1A1_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_2_2'): 'a2A1_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_3_0'): 'a0A12_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_3_1'): 'a1A12_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_3_2'): 'a2A12_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_4_0'): 'a0V_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_4_1'): 'a1V_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_4_2'): 'a2V_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_5_0'): 'a0T1_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_5_1'): 'a1T1_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_5_2'): 'a2T1_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_6_0'): 'a0T2_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_6_1'): 'a1T2_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_6_2'): 'a2T2_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_7_0'): 'a0T23_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_7_1'): 'a1T23_BKstar_BSZ_LCSRonly',
    ('B_Ks', '2_7_2'): 'a2T23_BKstar_BSZ_LCSRonly',
    ('B_Ks', '3_1_0'): 'a0A0_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_1_1'): 'a1A0_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_1_2'): 'a2A0_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_2_0'): 'a0A1_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_2_1'): 'a1A1_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_2_2'): 'a2A1_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_3_0'): 'a0A12_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_3_1'): 'a1A12_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_3_2'): 'a2A12_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_4_0'): 'a0V_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_4_1'): 'a1V_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_4_2'): 'a2V_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_5_0'): 'a0T1_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_5_1'): 'a1T1_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_5_2'): 'a2T1_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_6_0'): 'a0T2_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_6_1'): 'a1T2_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_6_2'): 'a2T2_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_7_0'): 'a0T23_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_7_1'): 'a1T23_BKstar_GRvDV_BSZ',
    ('B_Ks', '3_7_2'): 'a2T23_BKstar_GRvDV_BSZ',
    ('B_Ks', '4_1_0'): 'a0A0_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_1_1'): 'a1A0_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_1_2'): 'a2A0_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_2_0'): 'a0A1_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_2_1'): 'a1A1_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_2_2'): 'a2A1_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_3_0'): 'a0A12_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_3_1'): 'a1A12_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_3_2'): 'a2A12_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_4_0'): 'a0V_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_4_1'): 'a1V_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_4_2'): 'a2V_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_5_0'): 'a0T1_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_5_1'): 'a1T1_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_5_2'): 'a2T1_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_6_0'): 'a0T2_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_6_1'): 'a1T2_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_6_2'): 'a2T2_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_7_0'): 'a0T23_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_7_1'): 'a1T23_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '4_7_2'): 'a2T23_BKstar_GKvD_LCSR_Lattice',
    ('B_Ks', '5_1_0'): 'a0A0_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_1_1'): 'a1A0_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_1_2'): 'a2A0_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_2_0'): 'a0A1_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_2_1'): 'a1A1_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_2_2'): 'a2A1_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_3_0'): 'a0A12_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_3_1'): 'a1A12_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_3_2'): 'a2A12_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_4_0'): 'a0V_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_4_1'): 'a1V_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_4_2'): 'a2V_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_5_0'): 'a0T1_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_5_1'): 'a1T1_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_5_2'): 'a2T1_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_6_0'): 'a0T2_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_6_1'): 'a1T2_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_6_2'): 'a2T2_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_7_0'): 'a0T23_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_7_1'): 'a1T23_BKstar_GKvD_LCSRonly',
    ('B_Ks', '5_7_2'): 'a2T23_BKstar_GKvD_LCSRonly',
    ('B_Ks', '6'): 'systErr_BKstar_HLMW',
    ('B_Ks', '6_1_0'): 'a0A0_BKstar_HLMW',
    ('B_Ks', '6_1_1'): 'a1A0_BKstar_HLMW',
    ('B_Ks', '6_2_0'): 'a0A1_BKstar_HLMW',
    ('B_Ks', '6_2_1'): 'a1A1_BKstar_HLMW',
    ('B_Ks', '6_3_0'): 'a0A12_BKstar_HLMW',
    ('B_Ks', '6_3_1'): 'a1A12_BKstar_HLMW',
    ('B_Ks', '6_4_0'): 'a0V_BKstar_HLMW',
    ('B_Ks', '6_4_1'): 'a1V_BKstar_HLMW',
    ('B_Ks', '6_5_0'): 'a0T1_BKstar_HLMW',
    ('B_Ks', '6_5_1'): 'a1T1_BKstar_HLMW',
    ('B_Ks', '6_6_0'): 'a0T2_BKstar_HLMW',
    ('B_Ks', '6_6_1'): 'a1T2_BKstar_HLMW',
    ('B_Ks', '6_7_0'): 'a0T23_BKstar_HLMW',
    ('B_Ks', '6_7_1'): 'a1T23_BKstar_HLMW',
    ('B_Ks', '19_1_1_0'): 'real_alpha_perp0',
    ('B_Ks', '19_1_1_1'): 'real_alpha_perp1',
    ('B_Ks', '19_1_1_2'): 'real_alpha_perp2',
    ('B_Ks', '19_1_2_0'): 'real_alpha_par0',
    ('B_Ks', '19_1_2_1'): 'real_alpha_par1',
    ('B_Ks', '19_1_2_2'): 'real_alpha_par2',
    ('B_Ks', '19_1_3_0'): 'real_alpha_zero0',
    ('B_Ks', '19_1_3_1'): 'real_alpha_zero1',
    ('B_Ks', '19_2_1_0'): 'imag_alpha_perp0',
    ('B_Ks', '19_2_1_1'): 'imag_alpha_perp1',
    ('B_Ks', '19_2_1_2'): 'imag_alpha_perp2',
    ('B_Ks', '19_2_2_0'): 'imag_alpha_par0',
    ('B_Ks', '19_2_2_1'): 'imag_alpha_par1',
    ('B_Ks', '19_2_2_2'): 'imag_alpha_par2',
    ('B_Ks', '19_2_3_0'): 'imag_alpha_zero0',
    ('B_Ks', '19_2_3_1'): 'imag_alpha_zero1',
    ('B_Ks', '20_1'): 'DeltaC9_M1_q2bar',
    ('B_Ks', '20_2'): 'DeltaC9_M2_q2bar',
    ('B_Ks', '20_3'): 'DeltaC9_M3_q2bar',
    ('B_Ks', '21_1_1'): 'r1_M1',
    ('B_Ks', '21_1_2'): 'r1_M2',
    ('B_Ks', '21_1_3'): 'r1_M3',
    ('B_Ks', '21_2_1'): 'r2_M1',
    ('B_Ks', '21_2_2'): 'r2_M2',
    ('B_Ks', '21_2_3'): 'r2_M3',
    ('B_K', '1_1_0'): 'a0fp_BK_AS_LCSR_Lattice',
    ('B_K', '1_1_1'): 'a1fp_BK_AS_LCSR_Lattice',
    ('B_K', '1_1_2'): 'a2fp_BK_AS_LCSR_Lattice',
    ('B_K', '1_2_0'): 'a0f0_BK_AS_LCSR_Lattice',
    ('B_K', '1_2_1'): 'a1f0_BK_AS_LCSR_Lattice',
    ('B_K', '1_2_2'): 'a2f0_BK_AS_LCSR_Lattice',
    ('B_K', '1_2_3'): 'a3f0_BK_AS_LCSR_Lattice',
    ('B_K', '1_3_0'): 'a0fT_BK_AS_LCSR_Lattice',
    ('B_K', '1_3_1'): 'a1fT_BK_AS_LCSR_Lattice',
    ('B_K', '1_3_2'): 'a2fT_BK_AS_LCSR_Lattice',
    ('B_K', '2_1_0'): 'a0fp_BK_GRvDV_BSZ',
    ('B_K', '2_1_1'): 'a1fp_BK_GRvDV_BSZ',
    ('B_K', '2_1_2'): 'a2fp_BK_GRvDV_BSZ',
    ('B_K', '2_2_1'): 'a1f0_BK_GRvDV_BSZ',
    ('B_K', '2_2_2'): 'a2f0_BK_GRvDV_BSZ',
    ('B_K', '2_3_0'): 'a0fT_BK_GRvDV_BSZ',
    ('B_K', '2_3_1'): 'a1fT_BK_GRvDV_BSZ',
    ('B_K', '2_3_2'): 'a2fT_BK_GRvDV_BSZ',
    ('B_K', '3_1_0'): 'a0fp_BK_GKvD_LCSR_Lattice',
    ('B_K', '3_1_1'): 'a1fp_BK_GKvD_LCSR_Lattice',
    ('B_K', '3_1_2'): 'a2fp_BK_GKvD_LCSR_Lattice',
    ('B_K', '3_2_1'): 'a1f0_BK_GKvD_LCSR_Lattice',
    ('B_K', '3_2_2'): 'a2f0_BK_GKvD_LCSR_Lattice',
    ('B_K', '3_3_0'): 'a0fT_BK_GKvD_LCSR_Lattice',
    ('B_K', '3_3_1'): 'a1fT_BK_GKvD_LCSR_Lattice',
    ('B_K', '3_3_2'): 'a2fT_BK_GKvD_LCSR_Lattice',
    ('B_K', '4_1_0'): 'a0fp_BK_GKvD_LCSRonly',
    ('B_K', '4_1_1'): 'a1fp_BK_GKvD_LCSRonly',
    ('B_K', '4_1_2'): 'a2fp_BK_GKvD_LCSRonly',
    ('B_K', '4_2_1'): 'a1f0_BK_GKvD_LCSRonly',
    ('B_K', '4_2_2'): 'a2f0_BK_GKvD_LCSRonly',
    ('B_K', '4_3_0'): 'a0fT_BK_GKvD_LCSRonly',
    ('B_K', '4_3_1'): 'a1fT_BK_GKvD_LCSRonly',
    ('B_K', '4_3_2'): 'a2fT_BK_GKvD_LCSRonly',
    ('B_K', '5_1_0'): 'a0fp_BK_FLAG24',
    ('B_K', '5_1_1'): 'a1fp_BK_FLAG24',
    ('B_K', '5_1_2'): 'a2fp_BK_FLAG24',
    ('B_K', '5_2_0'): 'a0f0_BK_FLAG24',
    ('B_K', '5_2_1'): 'a1f0_BK_FLAG24',
    ('B_K', '5_3_0'): 'a0fT_BK_FLAG24',
    ('B_K', '5_3_1'): 'a1fT_BK_FLAG24',
    ('B_K', '5_3_2'): 'a2fT_BK_FLAG24',
    ('B_K', '6_1_1'): 'a1fp_BK_HPQCD22',
    ('B_K', '6_1_2'): 'a2fp_BK_HPQCD22',
    ('B_K', '6_2_0'): 'a0f0_BK_HPQCD22',
    ('B_K', '6_2_1'): 'a1f0_BK_HPQCD22',
    ('B_K', '6_2_2'): 'a2f0_BK_HPQCD22',
    ('B_K', '6_3_0'): 'a0fT_BK_HPQCD22',
    ('B_K', '6_3_1'): 'a1fT_BK_HPQCD22',
    ('B_K', '6_3_2'): 'a2fT_BK_HPQCD22',
    ('B_K', '8_1'): 'a1K',
    ('B_K', '8_2'): 'a2K',
    ('FLIFE', '531'): 'life_Bs',
    ('FCONST', '531_1'): 'f_Bs',
    ('B_ll', '1'): 'ys_Bs',
    ('FCONST', '333_1'): 'f_phi_par',
    ('FCONST', '333_2'): 'f_phi_perp',
    ('B_phi', '7_1'): 'a1phi_perp',
    ('B_phi', '7_2'): 'a2phi_perp',
    ('B_phi', '8_1'): 'a1phi_par',
    ('B_phi', '8_2'): 'a2phi_par',
    ('B_phi', '1_1_0'): 'a0A0_Bsphi',
    ('B_phi', '1_1_1'): 'a1A0_Bsphi',
    ('B_phi', '1_1_2'): 'a2A0_Bsphi',
    ('B_phi', '1_2_0'): 'a0A1_Bsphi',
    ('B_phi', '1_2_1'): 'a1A1_Bsphi',
    ('B_phi', '1_2_2'): 'a2A1_Bsphi',
    ('B_phi', '1_3_0'): 'a0A12_Bsphi',
    ('B_phi', '1_3_1'): 'a1A12_Bsphi',
    ('B_phi', '1_3_2'): 'a2A12_Bsphi',
    ('B_phi', '1_4_0'): 'a0V_Bsphi',
    ('B_phi', '1_4_1'): 'a1V_Bsphi',
    ('B_phi', '1_4_2'): 'a2V_Bsphi',
    ('B_phi', '1_5_0'): 'a0T1_Bsphi',
    ('B_phi', '1_5_1'): 'a1T1_Bsphi',
    ('B_phi', '1_5_2'): 'a2T1_Bsphi',
    ('B_phi', '1_6_0'): 'a0T2_Bsphi',
    ('B_phi', '1_6_1'): 'a1T2_Bsphi',
    ('B_phi', '1_6_2'): 'a2T2_Bsphi',
    ('B_phi', '1_7_0'): 'a0T23_Bsphi',
    ('B_phi', '1_7_1'): 'a1T23_Bsphi',
    ('B_phi', '1_7_2'): 'a2T23_Bsphi',
    ('B_phi', '2_1_0'): 'a0A0_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_1_1'): 'a1A0_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_1_2'): 'a2A0_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_2_0'): 'a0A1_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_2_1'): 'a1A1_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_2_2'): 'a2A1_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_3_0'): 'a0A12_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_3_1'): 'a1A12_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_3_2'): 'a2A12_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_4_0'): 'a0V_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_4_1'): 'a1V_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_4_2'): 'a2V_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_5_0'): 'a0T1_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_5_1'): 'a1T1_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_5_2'): 'a2T1_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_6_0'): 'a0T2_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_6_1'): 'a1T2_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_6_2'): 'a2T2_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_7_0'): 'a0T23_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_7_1'): 'a1T23_Bsphi_BSZ_LCSRonly',
    ('B_phi', '2_7_2'): 'a2T23_Bsphi_BSZ_LCSRonly',
    ('B_phi', '3_1_0'): 'a0A0_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_1_1'): 'a1A0_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_1_2'): 'a2A0_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_2_0'): 'a0A1_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_2_1'): 'a1A1_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_2_2'): 'a2A1_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_3_0'): 'a0A12_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_3_1'): 'a1A12_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_3_2'): 'a2A12_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_4_0'): 'a0V_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_4_1'): 'a1V_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_4_2'): 'a2V_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_5_0'): 'a0T1_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_5_1'): 'a1T1_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_5_2'): 'a2T1_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_6_0'): 'a0T2_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_6_1'): 'a1T2_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_6_2'): 'a2T2_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_7_0'): 'a0T23_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_7_1'): 'a1T23_Bsphi_GRvDV_BSZ',
    ('B_phi', '3_7_2'): 'a2T23_Bsphi_GRvDV_BSZ',
    ('B_phi', '6'): 'systErr_Bsphi_HLMW',
    ('B_phi', '6_1_0'): 'a0A0_Bsphi_HLMW',
    ('B_phi', '6_1_1'): 'a1A0_Bsphi_HLMW',
    ('B_phi', '6_2_0'): 'a0A1_Bsphi_HLMW',
    ('B_phi', '6_2_1'): 'a1A1_Bsphi_HLMW',
    ('B_phi', '6_3_0'): 'a0A12_Bsphi_HLMW',
    ('B_phi', '6_3_1'): 'a1A12_Bsphi_HLMW',
    ('B_phi', '6_4_0'): 'a0V_Bsphi_HLMW',
    ('B_phi', '6_4_1'): 'a1V_Bsphi_HLMW',
    ('B_phi', '6_5_0'): 'a0T1_Bsphi_HLMW',
    ('B_phi', '6_5_1'): 'a1T1_Bsphi_HLMW',
    ('B_phi', '6_6_0'): 'a0T2_Bsphi_HLMW',
    ('B_phi', '6_6_1'): 'a1T2_Bsphi_HLMW',
    ('B_phi', '6_7_0'): 'a0T23_Bsphi_HLMW',
    ('B_phi', '6_7_1'): 'a1T23_Bsphi_HLMW',
    ('FLIFE', '5122'): 'life_Lb',
    ('Lb_L', '8'): 'alphaL_LbLll',
    ('Lb_L', '1_1_0'): 'a0_HO_fperp_LbLll',
    ('Lb_L', '1_1_1'): 'a1_HO_fperp_LbLll',
    ('Lb_L', '1_1_2'): 'a2_HO_fperp_LbLll',
    ('Lb_L', '1_2_0'): 'a0_HO_hperp_LbLll',
    ('Lb_L', '1_2_1'): 'a1_HO_hperp_LbLll',
    ('Lb_L', '1_2_2'): 'a2_HO_hperp_LbLll',
    ('Lb_L', '1_3_1'): 'a1_HO_gperp_LbLll',
    ('Lb_L', '1_3_2'): 'a2_HO_gperp_LbLll',
    ('Lb_L', '1_4_0'): 'a0_HO_fplus_LbLll',
    ('Lb_L', '1_4_1'): 'a1_HO_fplus_LbLll',
    ('Lb_L', '1_4_2'): 'a2_HO_fplus_LbLll',
    ('Lb_L', '1_5_0'): 'a0_HO_hplus_LbLll',
    ('Lb_L', '1_5_1'): 'a1_HO_hplus_LbLll',
    ('Lb_L', '1_5_2'): 'a2_HO_hplus_LbLll',
    ('Lb_L', '1_6_1'): 'a1_HO_gplus_LbLll',
    ('Lb_L', '1_6_2'): 'a2_HO_gplus_LbLll',
    ('Lb_L', '1_7_1'): 'a1_HO_htildeperp_LbLll',
    ('Lb_L', '1_7_2'): 'a2_HO_htildeperp_LbLll',
    ('Lb_L', '1_8_1'): 'a1_HO_htildeplus_LbLll',
    ('Lb_L', '1_8_2'): 'a2_HO_htildeplus_LbLll',
    ('K_ll', '1'): 'BR_KLgammagamma_exp',
    ('K_ll', '2'): 'BR_KSgammagamma_exp',
    ('B_Ks', '18_1_1'): 'BtoKstarlow_ALperp_err_noq2',
    ('B_Ks', '18_1_2'): 'BtoKstarlow_ARperp_err_noq2',
    ('B_Ks', '18_1_3'): 'BtoKstarlow_ALpar_err_noq2',
    ('B_Ks', '18_1_4'): 'BtoKstarlow_ARpar_err_noq2',
    ('B_Ks', '18_1_5'): 'BtoKstarlow_AL0_err_noq2',
    ('B_Ks', '18_1_6'): 'BtoKstarlow_AR0_err_noq2',
    ('B_Ks', '18_1_7'): 'BtoKstarlow_At_err_noq2',
    ('B_Ks', '18_1_8'): 'BtoKstarlow_AS_err_noq2',
    ('B_Ks', '18_2_1'): 'BtoKstarlow_ALperp_err_q2',
    ('B_Ks', '18_2_2'): 'BtoKstarlow_ARperp_err_q2',
    ('B_Ks', '18_2_3'): 'BtoKstarlow_ALpar_err_q2',
    ('B_Ks', '18_2_4'): 'BtoKstarlow_ARpar_err_q2',
    ('B_Ks', '18_2_5'): 'BtoKstarlow_AL0_err_q2',
    ('B_Ks', '18_2_6'): 'BtoKstarlow_AR0_err_q2',
    ('B_Ks', '18_2_7'): 'BtoKstarlow_At_err_q2',
    ('B_Ks', '18_2_8'): 'BtoKstarlow_AS_err_q2',
    ('B_Ks', '18_3_1'): 'BtoKstarhigh_ALperp_err',
    ('B_Ks', '18_3_2'): 'BtoKstarhigh_ARperp_err',
    ('B_Ks', '18_3_3'): 'BtoKstarhigh_ALpar_err',
    ('B_Ks', '18_3_4'): 'BtoKstarhigh_ARpar_err',
    ('B_Ks', '18_3_5'): 'BtoKstarhigh_AL0_err',
    ('B_Ks', '18_3_6'): 'BtoKstarhigh_AR0_err',
    ('B_Ks', '18_3_7'): 'BtoKstarhigh_At_err',
    ('B_Ks', '18_3_8'): 'BtoKstarhigh_AS_err',
    ('B_K', '18_1_1'): 'BtoKlow_FV_err_noq2',
    ('B_K', '18_1_2'): 'BtoKlow_FA_err_noq2',
    ('B_K', '18_1_3'): 'BtoKlow_FS_err_noq2',
    ('B_K', '18_1_4'): 'BtoKlow_FP_err_noq2',
    ('B_K', '18_2_1'): 'BtoKlow_FV_err_q2',
    ('B_K', '18_2_2'): 'BtoKlow_FA_err_q2',
    ('B_K', '18_2_3'): 'BtoKlow_FS_err_q2',
    ('B_K', '18_2_4'): 'BtoKlow_FP_err_q2',
    ('B_K', '18_3_1'): 'BtoKhigh_FV_err',
    ('B_K', '18_3_2'): 'BtoKhigh_FA_err',
    ('B_K', '18_3_3'): 'BtoKhigh_FS_err',
    ('B_K', '18_3_4'): 'BtoKhigh_FP_err',
    ('B_phi', '18_1_1'): 'Bstophilow_ALperp_err_noq2',
    ('B_phi', '18_1_2'): 'Bstophilow_ARperp_err_noq2',
    ('B_phi', '18_1_3'): 'Bstophilow_ALpar_err_noq2',
    ('B_phi', '18_1_4'): 'Bstophilow_ARpar_err_noq2',
    ('B_phi', '18_1_5'): 'Bstophilow_AL0_err_noq2',
    ('B_phi', '18_1_6'): 'Bstophilow_AR0_err_noq2',
    ('B_phi', '18_1_7'): 'Bstophilow_At_err_noq2',
    ('B_phi', '18_1_8'): 'Bstophilow_AS_err_noq2',
    ('B_phi', '18_2_1'): 'Bstophilow_ALperp_err_q2',
    ('B_phi', '18_2_2'): 'Bstophilow_ARperp_err_q2',
    ('B_phi', '18_2_3'): 'Bstophilow_ALpar_err_q2',
    ('B_phi', '18_2_4'): 'Bstophilow_ARpar_err_q2',
    ('B_phi', '18_2_5'): 'Bstophilow_AL0_err_q2',
    ('B_phi', '18_2_6'): 'Bstophilow_AR0_err_q2',
    ('B_phi', '18_2_7'): 'Bstophilow_At_err_q2',
    ('B_phi', '18_2_8'): 'Bstophilow_AS_err_q2',
    ('B_phi', '18_3_1'): 'Bstophihigh_ALperp_err',
    ('B_phi', '18_3_2'): 'Bstophihigh_ARperp_err',
    ('B_phi', '18_3_3'): 'Bstophihigh_ALpar_err',
    ('B_phi', '18_3_4'): 'Bstophihigh_ARpar_err',
    ('B_phi', '18_3_5'): 'Bstophihigh_AL0_err',
    ('B_phi', '18_3_6'): 'Bstophihigh_AR0_err',
    ('B_phi', '18_3_7'): 'Bstophihigh_At_err',
    ('B_phi', '18_3_8'): 'Bstophihigh_AS_err',
    ('FCONST', '511_1'): 'f_B',
    ('FCONST', '521_1'): 'f_B',
    ('FCONST', '311_1'): 'f_K',
    ('FCONST', '321_1'): 'f_K',
    ('B_K', '13'): 'lambda_Bp',
    ('B_Ks', '13'): 'lambda_Bp',
    ('B_phi', '13'): 'lambda_Bsp',
    ('Lb_L', '1_3_0'): 'a0_HO_gpp_LbLll',
    ('Lb_L', '1_6_0'): 'a0_HO_gpp_LbLll',
    ('Lb_L', '1_7_0'): 'a0_HO_htildepp_LbLll',
    ('Lb_L', '1_8_0'): 'a0_HO_htildepp_LbLll',
}

UNMATCHED_OR_NOT_EXPOSED_IN_NEW = [
    'log_mu_spec_lambda_h_mass_b',
    'log_muK_1GeV',
    'log_mu_W_mass_W',
    'log_mu_b_mass_b',
    'Aterm_mu_KLmumu',
    'Iterm_mu_KSmumu',
    'chi_gg_Mrho',
    'err_Pc_Xlambda_Kppipnunu',
    'deltaPcu_Kppipnunu',
    'KLpill_Ce_Cdir',
    'KLpill_Cmu_Cdir',
    'KLpill_Ce_Cint',
    'KLpill_Cmu_Cint',
    'KLpill_Ce_Cmix',
    'KLpill_Cmu_Cmix',
    'KLpill_Cmu_CPC',
    'KLpill_abs_aS',
    'KLpill_Ce_Sgg',
    'KLpill_Cmu_Sgg',
]


def get_latex_nuisance_name(block: str, pdg: str) -> str | None:
    """Return the LaTeX label for a new-format nuisance target."""
    return NEW_TO_LATEX_NUISANCE_MAP.get((block, pdg))


def get_legacy_nuisance_name(block: str, pdg: str) -> str | None:
    """Return the legacy nuisance id used to derive a new-format target."""
    return NEW_TO_LEGACY_NUISANCE_MAP.get((block, pdg))
