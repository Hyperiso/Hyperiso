#!/usr/bin/env python3
"""
Manual map between legacy nuisance IDs and the new (block, id) format.

Conventions
-----------
- Direct one-to-one mapping:
    LEGACY_TO_NEW_NUISANCE_MAP["alphas_MZ"] == ("SMINPUTS", "3")

- One legacy nuisance reused/split into several new entries:
    LEGACY_TO_NEW_MULTI_MAP["f_B"] == [("FCONST", "511_1"), ("FCONST", "521_1")]

- Several old nuisance names compressed onto a single new nuisance entry:
    COMPRESSED_OR_SHARED_MAP["BtoKstarlow_ALperp_err_noq2"] == ("B_Ks", "18_1")

- Approximate / best-effort matches (same physics role, but value or packaging changed):
    APPROXIMATE_MAP["mu_c_bsg"] == ("B_Xs", "7")

The goal is practical conversion help, not a proof that the statistical model is
strictly identical across formats.
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Union

NewTarget = Tuple[str, str]
NewTargets = List[NewTarget]


# ---------------------------------------------------------------------------
# 1) Strong one-to-one matches
# ---------------------------------------------------------------------------

LEGACY_TO_NEW_NUISANCE_MAP: Dict[str, NewTarget] = {
    # --- SM / CKM / masses / scales ---
    "alphas_MZ": ("SMINPUTS", "3"),
    "mass_b": ("SMINPUTS", "5"),
    "mass_c": ("MASS", "4"),
    "mass_s": ("MASS", "3"),
    "mass_top_pole": ("SMINPUTS", "6"),
    "mass_h0": ("MASS", "25"),
    "CKM_lambda": ("VCKMIN", "1"),
    "CKM_A": ("VCKMIN", "2"),
    "CKM_rhobar": ("VCKMIN", "3"),
    "CKM_etabar": ("VCKMIN", "4"),

    # --- inclusive B -> Xs ---
    "mu_c_bsg": ("B_Xs", "7"),
    "BR_BXclnu_exp": ("B_Xs", "2"),
    "mu_G2_bsg": ("B_Xs", "3"),
    "rho_D3_bsg": ("B_Xs", "4"),
    "rho_LS3_bsg": ("B_Xs", "5"),
    "bsgamma_rand": ("B_Xs", "10"),

    "BRBXstautau_lowq2_rand": ("B_Xsll", "3_15_0"),
    "BRBXstautau_highq2_rand": ("B_Xsll", "3_15_1"),
    "BRBXstautau_full_rand": ("B_Xsll", "3_15_2"),
    "BRBXsmumu_lowq2_rand": ("B_Xsll", "3_13_0"),
    "BRBXsmumu_highq2_rand": ("B_Xsll", "3_13_1"),
    "BRBXsmumu_full_rand": ("B_Xsll", "3_13_2"),
    "BRBXsee_lowq2_rand": ("B_Xsll", "3_11_0"),
    "BRBXsee_highq2_rand": ("B_Xsll", "3_11_1"),
    "BRBXsee_full_rand": ("B_Xsll", "3_11_2"),

    # --- B -> K* form factors / LCDAs / input constants ---
    "f_Kstar_par": ("FCONST", "323_1"),
    "f_Kstar_perp": ("FCONST", "323_2"),
    "a1perp": ("B_Ks", "7_1"),
    "a2perp": ("B_Ks", "7_2"),
    "a1par": ("B_Ks", "8_1"),
    "a2par": ("B_Ks", "8_2"),
    "T1_BKstar": ("B_Ks", "16"),

    "a0A0_BKstar": ("B_Ks", "1_1_0"),
    "a1A0_BKstar": ("B_Ks", "1_1_1"),
    "a2A0_BKstar": ("B_Ks", "1_1_2"),
    "a0A1_BKstar": ("B_Ks", "1_2_0"),
    "a1A1_BKstar": ("B_Ks", "1_2_1"),
    "a2A1_BKstar": ("B_Ks", "1_2_2"),
    "a0A12_BKstar": ("B_Ks", "1_3_0"),
    "a1A12_BKstar": ("B_Ks", "1_3_1"),
    "a2A12_BKstar": ("B_Ks", "1_3_2"),
    "a0V_BKstar": ("B_Ks", "1_4_0"),
    "a1V_BKstar": ("B_Ks", "1_4_1"),
    "a2V_BKstar": ("B_Ks", "1_4_2"),
    "a0T1_BKstar": ("B_Ks", "1_5_0"),
    "a1T1_BKstar": ("B_Ks", "1_5_1"),
    "a2T1_BKstar": ("B_Ks", "1_5_2"),
    "a0T2_BKstar": ("B_Ks", "1_6_0"),
    "a1T2_BKstar": ("B_Ks", "1_6_1"),
    "a2T2_BKstar": ("B_Ks", "1_6_2"),
    "a0T23_BKstar": ("B_Ks", "1_7_0"),
    "a1T23_BKstar": ("B_Ks", "1_7_1"),
    "a2T23_BKstar": ("B_Ks", "1_7_2"),

    "a0A0_BKstar_BSZ_LCSRonly": ("B_Ks", "2_1_0"),
    "a1A0_BKstar_BSZ_LCSRonly": ("B_Ks", "2_1_1"),
    "a2A0_BKstar_BSZ_LCSRonly": ("B_Ks", "2_1_2"),
    "a0A1_BKstar_BSZ_LCSRonly": ("B_Ks", "2_2_0"),
    "a1A1_BKstar_BSZ_LCSRonly": ("B_Ks", "2_2_1"),
    "a2A1_BKstar_BSZ_LCSRonly": ("B_Ks", "2_2_2"),
    "a0A12_BKstar_BSZ_LCSRonly": ("B_Ks", "2_3_0"),
    "a1A12_BKstar_BSZ_LCSRonly": ("B_Ks", "2_3_1"),
    "a2A12_BKstar_BSZ_LCSRonly": ("B_Ks", "2_3_2"),
    "a0V_BKstar_BSZ_LCSRonly": ("B_Ks", "2_4_0"),
    "a1V_BKstar_BSZ_LCSRonly": ("B_Ks", "2_4_1"),
    "a2V_BKstar_BSZ_LCSRonly": ("B_Ks", "2_4_2"),
    "a0T1_BKstar_BSZ_LCSRonly": ("B_Ks", "2_5_0"),
    "a1T1_BKstar_BSZ_LCSRonly": ("B_Ks", "2_5_1"),
    "a2T1_BKstar_BSZ_LCSRonly": ("B_Ks", "2_5_2"),
    "a0T2_BKstar_BSZ_LCSRonly": ("B_Ks", "2_6_0"),
    "a1T2_BKstar_BSZ_LCSRonly": ("B_Ks", "2_6_1"),
    "a2T2_BKstar_BSZ_LCSRonly": ("B_Ks", "2_6_2"),
    "a0T23_BKstar_BSZ_LCSRonly": ("B_Ks", "2_7_0"),
    "a1T23_BKstar_BSZ_LCSRonly": ("B_Ks", "2_7_1"),
    "a2T23_BKstar_BSZ_LCSRonly": ("B_Ks", "2_7_2"),

    "a0A0_BKstar_GRvDV_BSZ": ("B_Ks", "3_1_0"),
    "a1A0_BKstar_GRvDV_BSZ": ("B_Ks", "3_1_1"),
    "a2A0_BKstar_GRvDV_BSZ": ("B_Ks", "3_1_2"),
    "a0A1_BKstar_GRvDV_BSZ": ("B_Ks", "3_2_0"),
    "a1A1_BKstar_GRvDV_BSZ": ("B_Ks", "3_2_1"),
    "a2A1_BKstar_GRvDV_BSZ": ("B_Ks", "3_2_2"),
    "a0A12_BKstar_GRvDV_BSZ": ("B_Ks", "3_3_0"),
    "a1A12_BKstar_GRvDV_BSZ": ("B_Ks", "3_3_1"),
    "a2A12_BKstar_GRvDV_BSZ": ("B_Ks", "3_3_2"),
    "a0V_BKstar_GRvDV_BSZ": ("B_Ks", "3_4_0"),
    "a1V_BKstar_GRvDV_BSZ": ("B_Ks", "3_4_1"),
    "a2V_BKstar_GRvDV_BSZ": ("B_Ks", "3_4_2"),
    "a0T1_BKstar_GRvDV_BSZ": ("B_Ks", "3_5_0"),
    "a1T1_BKstar_GRvDV_BSZ": ("B_Ks", "3_5_1"),
    "a2T1_BKstar_GRvDV_BSZ": ("B_Ks", "3_5_2"),
    "a0T2_BKstar_GRvDV_BSZ": ("B_Ks", "3_6_0"),
    "a1T2_BKstar_GRvDV_BSZ": ("B_Ks", "3_6_1"),
    "a2T2_BKstar_GRvDV_BSZ": ("B_Ks", "3_6_2"),
    "a0T23_BKstar_GRvDV_BSZ": ("B_Ks", "3_7_0"),
    "a1T23_BKstar_GRvDV_BSZ": ("B_Ks", "3_7_1"),
    "a2T23_BKstar_GRvDV_BSZ": ("B_Ks", "3_7_2"),

    "a0A0_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_1_0"),
    "a1A0_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_1_1"),
    "a2A0_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_1_2"),
    "a0A1_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_2_0"),
    "a1A1_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_2_1"),
    "a2A1_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_2_2"),
    "a0A12_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_3_0"),
    "a1A12_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_3_1"),
    "a2A12_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_3_2"),
    "a0V_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_4_0"),
    "a1V_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_4_1"),
    "a2V_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_4_2"),
    "a0T1_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_5_0"),
    "a1T1_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_5_1"),
    "a2T1_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_5_2"),
    "a0T2_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_6_0"),
    "a1T2_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_6_1"),
    "a2T2_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_6_2"),
    "a0T23_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_7_0"),
    "a1T23_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_7_1"),
    "a2T23_BKstar_GKvD_LCSR_Lattice": ("B_Ks", "4_7_2"),

    "a0A0_BKstar_GKvD_LCSRonly": ("B_Ks", "5_1_0"),
    "a1A0_BKstar_GKvD_LCSRonly": ("B_Ks", "5_1_1"),
    "a2A0_BKstar_GKvD_LCSRonly": ("B_Ks", "5_1_2"),
    "a0A1_BKstar_GKvD_LCSRonly": ("B_Ks", "5_2_0"),
    "a1A1_BKstar_GKvD_LCSRonly": ("B_Ks", "5_2_1"),
    "a2A1_BKstar_GKvD_LCSRonly": ("B_Ks", "5_2_2"),
    "a0A12_BKstar_GKvD_LCSRonly": ("B_Ks", "5_3_0"),
    "a1A12_BKstar_GKvD_LCSRonly": ("B_Ks", "5_3_1"),
    "a2A12_BKstar_GKvD_LCSRonly": ("B_Ks", "5_3_2"),
    "a0V_BKstar_GKvD_LCSRonly": ("B_Ks", "5_4_0"),
    "a1V_BKstar_GKvD_LCSRonly": ("B_Ks", "5_4_1"),
    "a2V_BKstar_GKvD_LCSRonly": ("B_Ks", "5_4_2"),
    "a0T1_BKstar_GKvD_LCSRonly": ("B_Ks", "5_5_0"),
    "a1T1_BKstar_GKvD_LCSRonly": ("B_Ks", "5_5_1"),
    "a2T1_BKstar_GKvD_LCSRonly": ("B_Ks", "5_5_2"),
    "a0T2_BKstar_GKvD_LCSRonly": ("B_Ks", "5_6_0"),
    "a1T2_BKstar_GKvD_LCSRonly": ("B_Ks", "5_6_1"),
    "a2T2_BKstar_GKvD_LCSRonly": ("B_Ks", "5_6_2"),
    "a0T23_BKstar_GKvD_LCSRonly": ("B_Ks", "5_7_0"),
    "a1T23_BKstar_GKvD_LCSRonly": ("B_Ks", "5_7_1"),
    "a2T23_BKstar_GKvD_LCSRonly": ("B_Ks", "5_7_2"),

    "systErr_BKstar_HLMW" : ("B_Ks", "6"),
    "a0A0_BKstar_HLMW": ("B_Ks", "6_1_0"),
    "a1A0_BKstar_HLMW": ("B_Ks", "6_1_1"),
    "a0A1_BKstar_HLMW": ("B_Ks", "6_2_0"),
    "a1A1_BKstar_HLMW": ("B_Ks", "6_2_1"),
    "a0A12_BKstar_HLMW": ("B_Ks", "6_3_0"),
    "a1A12_BKstar_HLMW": ("B_Ks", "6_3_1"),
    "a0V_BKstar_HLMW": ("B_Ks", "6_4_0"),
    "a1V_BKstar_HLMW": ("B_Ks", "6_4_1"),
    "a0T1_BKstar_HLMW": ("B_Ks", "6_5_0"),
    "a1T1_BKstar_HLMW": ("B_Ks", "6_5_1"),
    "a0T2_BKstar_HLMW": ("B_Ks", "6_6_0"),
    "a1T2_BKstar_HLMW": ("B_Ks", "6_6_1"),
    "a0T23_BKstar_HLMW": ("B_Ks", "6_7_0"),
    "a1T23_BKstar_HLMW": ("B_Ks", "6_7_1"),

    "real_alpha_perp0" : ("B_Ks", "19_1_1_0"),
    "real_alpha_perp1" : ("B_Ks", "19_1_1_1"),
    "real_alpha_perp2" : ("B_Ks", "19_1_1_2"),
    "real_alpha_par0" : ("B_Ks", "19_1_2_0"),
    "real_alpha_par1" : ("B_Ks", "19_1_2_1"),
    "real_alpha_par2" : ("B_Ks", "19_1_2_2"),
    "real_alpha_zero0" : ("B_Ks", "19_1_3_0"),
    "real_alpha_zero1" : ("B_Ks", "19_1_3_1"),
    "imag_alpha_perp0" : ("B_Ks", "19_2_1_0"),
    "imag_alpha_perp1" : ("B_Ks", "19_2_1_1"),
    "imag_alpha_perp2" : ("B_Ks", "19_2_1_2"),
    "imag_alpha_par0" : ("B_Ks", "19_2_2_0"),
    "imag_alpha_par1" : ("B_Ks", "19_2_2_1"),
    "imag_alpha_par2" : ("B_Ks", "19_2_2_2"),
    "imag_alpha_zero0" : ("B_Ks", "19_2_3_0"),
    "imag_alpha_zero1" : ("B_Ks", "19_2_3_1"),

    "DeltaC9_M1_q2bar" : ("B_Ks", "20_1"),
    "DeltaC9_M2_q2bar" : ("B_Ks", "20_2"),
    "DeltaC9_M3_q2bar" : ("B_Ks", "20_3"),
    "r1_M1" : ("B_Ks", "21_1_1"),
    "r1_M2" : ("B_Ks", "21_1_2"),
    "r1_M3" : ("B_Ks", "21_1_3"),
    "r2_M1" : ("B_Ks", "21_2_1"),
    "r2_M2" : ("B_Ks", "21_2_2"),
    "r2_M3" : ("B_Ks", "21_2_3"),


    # --- B -> K form factors / LCDAs ---
    "a0fp_BK_AS_LCSR_Lattice": ("B_K", "1_1_0"),
    "a1fp_BK_AS_LCSR_Lattice": ("B_K", "1_1_1"),
    "a2fp_BK_AS_LCSR_Lattice": ("B_K", "1_1_2"),
    "a0f0_BK_AS_LCSR_Lattice": ("B_K", "1_2_0"),
    "a1f0_BK_AS_LCSR_Lattice": ("B_K", "1_2_1"),
    "a2f0_BK_AS_LCSR_Lattice": ("B_K", "1_2_2"),
    "a3f0_BK_AS_LCSR_Lattice": ("B_K", "1_2_3"),
    "a0fT_BK_AS_LCSR_Lattice": ("B_K", "1_3_0"),
    "a1fT_BK_AS_LCSR_Lattice": ("B_K", "1_3_1"),
    "a2fT_BK_AS_LCSR_Lattice": ("B_K", "1_3_2"),

    "a0fp_BK_GRvDV_BSZ": ("B_K", "2_1_0"),
    "a1fp_BK_GRvDV_BSZ": ("B_K", "2_1_1"),
    "a2fp_BK_GRvDV_BSZ": ("B_K", "2_1_2"),
    "a1f0_BK_GRvDV_BSZ": ("B_K", "2_2_1"),
    "a2f0_BK_GRvDV_BSZ": ("B_K", "2_2_2"),
    "a0fT_BK_GRvDV_BSZ": ("B_K", "2_3_0"),
    "a1fT_BK_GRvDV_BSZ": ("B_K", "2_3_1"),
    "a2fT_BK_GRvDV_BSZ": ("B_K", "2_3_2"),

    "a0fp_BK_GKvD_LCSR_Lattice": ("B_K", "3_1_0"),
    "a1fp_BK_GKvD_LCSR_Lattice": ("B_K", "3_1_1"),
    "a2fp_BK_GKvD_LCSR_Lattice": ("B_K", "3_1_2"),
    "a1f0_BK_GKvD_LCSR_Lattice": ("B_K", "3_2_1"),
    "a2f0_BK_GKvD_LCSR_Lattice": ("B_K", "3_2_2"),
    "a0fT_BK_GKvD_LCSR_Lattice": ("B_K", "3_3_0"),
    "a1fT_BK_GKvD_LCSR_Lattice": ("B_K", "3_3_1"),
    "a2fT_BK_GKvD_LCSR_Lattice": ("B_K", "3_3_2"),

    "a0fp_BK_GKvD_LCSRonly": ("B_K", "4_1_0"),
    "a1fp_BK_GKvD_LCSRonly": ("B_K", "4_1_1"),
    "a2fp_BK_GKvD_LCSRonly": ("B_K", "4_1_2"),
    "a1f0_BK_GKvD_LCSRonly": ("B_K", "4_2_1"),
    "a2f0_BK_GKvD_LCSRonly": ("B_K", "4_2_2"),
    "a0fT_BK_GKvD_LCSRonly": ("B_K", "4_3_0"),
    "a1fT_BK_GKvD_LCSRonly": ("B_K", "4_3_1"),
    "a2fT_BK_GKvD_LCSRonly": ("B_K", "4_3_2"),

    "a0fp_BK_FLAG24": ("B_K", "5_1_0"),
    "a1fp_BK_FLAG24": ("B_K", "5_1_1"),
    "a2fp_BK_FLAG24": ("B_K", "5_1_2"),
    "a0f0_BK_FLAG24": ("B_K", "5_2_0"),
    "a1f0_BK_FLAG24": ("B_K", "5_2_1"),
    "a0fT_BK_FLAG24": ("B_K", "5_3_0"),
    "a1fT_BK_FLAG24": ("B_K", "5_3_1"),
    "a2fT_BK_FLAG24": ("B_K", "5_3_2"),

    "a1fp_BK_HPQCD22": ("B_K", "6_1_1"),
    "a2fp_BK_HPQCD22": ("B_K", "6_1_2"),
    "a0f0_BK_HPQCD22": ("B_K", "6_2_0"),
    "a1f0_BK_HPQCD22": ("B_K", "6_2_1"),
    "a2f0_BK_HPQCD22": ("B_K", "6_2_2"),
    "a0fT_BK_HPQCD22": ("B_K", "6_3_0"),
    "a1fT_BK_HPQCD22": ("B_K", "6_3_1"),
    "a2fT_BK_HPQCD22": ("B_K", "6_3_2"),

    "a1K": ("B_K", "8_1"),
    "a2K": ("B_K", "8_2"),

    # --- Bs -> phi ---
    "life_Bs": ("FLIFE", "531"),
    "f_Bs": ("FCONST", "531_1"),
    "ys_Bs": ("B_ll", "1"),
    "f_phi_par": ("FCONST", "333_1"),
    "f_phi_perp": ("FCONST", "333_2"),
    "a1phi_perp": ("B_phi", "7_1"),
    "a2phi_perp": ("B_phi", "7_2"),
    "a1phi_par": ("B_phi", "8_1"),
    "a2phi_par": ("B_phi", "8_2"),

    "a0A0_Bsphi": ("B_phi", "1_1_0"),
    "a1A0_Bsphi": ("B_phi", "1_1_1"),
    "a2A0_Bsphi": ("B_phi", "1_1_2"),
    "a0A1_Bsphi": ("B_phi", "1_2_0"),
    "a1A1_Bsphi": ("B_phi", "1_2_1"),
    "a2A1_Bsphi": ("B_phi", "1_2_2"),
    "a0A12_Bsphi": ("B_phi", "1_3_0"),
    "a1A12_Bsphi": ("B_phi", "1_3_1"),
    "a2A12_Bsphi": ("B_phi", "1_3_2"),
    "a0V_Bsphi": ("B_phi", "1_4_0"),
    "a1V_Bsphi": ("B_phi", "1_4_1"),
    "a2V_Bsphi": ("B_phi", "1_4_2"),
    "a0T1_Bsphi": ("B_phi", "1_5_0"),
    "a1T1_Bsphi": ("B_phi", "1_5_1"),
    "a2T1_Bsphi": ("B_phi", "1_5_2"),
    "a0T2_Bsphi": ("B_phi", "1_6_0"),
    "a1T2_Bsphi": ("B_phi", "1_6_1"),
    "a2T2_Bsphi": ("B_phi", "1_6_2"),
    "a0T23_Bsphi": ("B_phi", "1_7_0"),
    "a1T23_Bsphi": ("B_phi", "1_7_1"),
    "a2T23_Bsphi": ("B_phi", "1_7_2"),

    "a0A0_Bsphi_BSZ_LCSRonly": ("B_phi", "2_1_0"),
    "a1A0_Bsphi_BSZ_LCSRonly": ("B_phi", "2_1_1"),
    "a2A0_Bsphi_BSZ_LCSRonly": ("B_phi", "2_1_2"),
    "a0A1_Bsphi_BSZ_LCSRonly": ("B_phi", "2_2_0"),
    "a1A1_Bsphi_BSZ_LCSRonly": ("B_phi", "2_2_1"),
    "a2A1_Bsphi_BSZ_LCSRonly": ("B_phi", "2_2_2"),
    "a0A12_Bsphi_BSZ_LCSRonly": ("B_phi", "2_3_0"),
    "a1A12_Bsphi_BSZ_LCSRonly": ("B_phi", "2_3_1"),
    "a2A12_Bsphi_BSZ_LCSRonly": ("B_phi", "2_3_2"),
    "a0V_Bsphi_BSZ_LCSRonly": ("B_phi", "2_4_0"),
    "a1V_Bsphi_BSZ_LCSRonly": ("B_phi", "2_4_1"),
    "a2V_Bsphi_BSZ_LCSRonly": ("B_phi", "2_4_2"),
    "a0T1_Bsphi_BSZ_LCSRonly": ("B_phi", "2_5_0"),
    "a1T1_Bsphi_BSZ_LCSRonly": ("B_phi", "2_5_1"),
    "a2T1_Bsphi_BSZ_LCSRonly": ("B_phi", "2_5_2"),
    "a0T2_Bsphi_BSZ_LCSRonly": ("B_phi", "2_6_0"),
    "a1T2_Bsphi_BSZ_LCSRonly": ("B_phi", "2_6_1"),
    "a2T2_Bsphi_BSZ_LCSRonly": ("B_phi", "2_6_2"),
    "a0T23_Bsphi_BSZ_LCSRonly": ("B_phi", "2_7_0"),
    "a1T23_Bsphi_BSZ_LCSRonly": ("B_phi", "2_7_1"),
    "a2T23_Bsphi_BSZ_LCSRonly": ("B_phi", "2_7_2"),

    "a0A0_Bsphi_GRvDV_BSZ": ("B_phi", "3_1_0"),
    "a1A0_Bsphi_GRvDV_BSZ": ("B_phi", "3_1_1"),
    "a2A0_Bsphi_GRvDV_BSZ": ("B_phi", "3_1_2"),
    "a0A1_Bsphi_GRvDV_BSZ": ("B_phi", "3_2_0"),
    "a1A1_Bsphi_GRvDV_BSZ": ("B_phi", "3_2_1"),
    "a2A1_Bsphi_GRvDV_BSZ": ("B_phi", "3_2_2"),
    "a0A12_Bsphi_GRvDV_BSZ": ("B_phi", "3_3_0"),
    "a1A12_Bsphi_GRvDV_BSZ": ("B_phi", "3_3_1"),
    "a2A12_Bsphi_GRvDV_BSZ": ("B_phi", "3_3_2"),
    "a0V_Bsphi_GRvDV_BSZ": ("B_phi", "3_4_0"),
    "a1V_Bsphi_GRvDV_BSZ": ("B_phi", "3_4_1"),
    "a2V_Bsphi_GRvDV_BSZ": ("B_phi", "3_4_2"),
    "a0T1_Bsphi_GRvDV_BSZ": ("B_phi", "3_5_0"),
    "a1T1_Bsphi_GRvDV_BSZ": ("B_phi", "3_5_1"),
    "a2T1_Bsphi_GRvDV_BSZ": ("B_phi", "3_5_2"),
    "a0T2_Bsphi_GRvDV_BSZ": ("B_phi", "3_6_0"),
    "a1T2_Bsphi_GRvDV_BSZ": ("B_phi", "3_6_1"),
    "a2T2_Bsphi_GRvDV_BSZ": ("B_phi", "3_6_2"),
    "a0T23_Bsphi_GRvDV_BSZ": ("B_phi", "3_7_0"),
    "a1T23_Bsphi_GRvDV_BSZ": ("B_phi", "3_7_1"),
    "a2T23_Bsphi_GRvDV_BSZ": ("B_phi", "3_7_2"),

    "systErr_Bsphi_HLMW" : ("B_phi", "6"),
    "a0A0_Bsphi_HLMW": ("B_phi", "6_1_0"),
    "a1A0_Bsphi_HLMW": ("B_phi", "6_1_1"),
    "a0A1_Bsphi_HLMW": ("B_phi", "6_2_0"),
    "a1A1_Bsphi_HLMW": ("B_phi", "6_2_1"),
    "a0A12_Bsphi_HLMW": ("B_phi", "6_3_0"),
    "a1A12_Bsphi_HLMW": ("B_phi", "6_3_1"),
    "a0V_Bsphi_HLMW": ("B_phi", "6_4_0"),
    "a1V_Bsphi_HLMW": ("B_phi", "6_4_1"),
    "a0T1_Bsphi_HLMW": ("B_phi", "6_5_0"),
    "a1T1_Bsphi_HLMW": ("B_phi", "6_5_1"),
    "a0T2_Bsphi_HLMW": ("B_phi", "6_6_0"),
    "a1T2_Bsphi_HLMW": ("B_phi", "6_6_1"),
    "a0T23_Bsphi_HLMW": ("B_phi", "6_7_0"),
    "a1T23_Bsphi_HLMW": ("B_phi", "6_7_1"),

    # --- Lambda_b -> Lambda ---
    "life_Lb": ("FLIFE", "5122"),
    "alphaL_LbLll": ("Lb_L", "8"),
    "a0_HO_fperp_LbLll": ("Lb_L", "1_1_0"),
    "a1_HO_fperp_LbLll": ("Lb_L", "1_1_1"),
    "a2_HO_fperp_LbLll": ("Lb_L", "1_1_2"),
    "a0_HO_hperp_LbLll": ("Lb_L", "1_2_0"),
    "a1_HO_hperp_LbLll": ("Lb_L", "1_2_1"),
    "a2_HO_hperp_LbLll": ("Lb_L", "1_2_2"),
    "a1_HO_gperp_LbLll": ("Lb_L", "1_3_1"),
    "a2_HO_gperp_LbLll": ("Lb_L", "1_3_2"),
    "a0_HO_fplus_LbLll": ("Lb_L", "1_4_0"),
    "a1_HO_fplus_LbLll": ("Lb_L", "1_4_1"),
    "a2_HO_fplus_LbLll": ("Lb_L", "1_4_2"),
    "a0_HO_hplus_LbLll": ("Lb_L", "1_5_0"),
    "a1_HO_hplus_LbLll": ("Lb_L", "1_5_1"),
    "a2_HO_hplus_LbLll": ("Lb_L", "1_5_2"),
    "a1_HO_gplus_LbLll": ("Lb_L", "1_6_1"),
    "a2_HO_gplus_LbLll": ("Lb_L", "1_6_2"),
    "a1_HO_htildeperp_LbLll": ("Lb_L", "1_7_1"),
    "a2_HO_htildeperp_LbLll": ("Lb_L", "1_7_2"),
    "a1_HO_htildeplus_LbLll": ("Lb_L", "1_8_1"),
    "a2_HO_htildeplus_LbLll": ("Lb_L", "1_8_2"),

    # --- Kaon / rare-K auxiliary inputs ---
    "BR_KLgammagamma_exp": ("K_ll", "1"),
    "BR_KSgammagamma_exp": ("K_ll", "2"),

    # --- B -> K* low-q2 non-factorizable nuisances ---
    "BtoKstarlow_ALperp_err_noq2": ("B_Ks", "18_1_1"),
    "BtoKstarlow_ARperp_err_noq2": ("B_Ks", "18_1_2"),
    "BtoKstarlow_ALpar_err_noq2": ("B_Ks", "18_1_3"),
    "BtoKstarlow_ARpar_err_noq2": ("B_Ks", "18_1_4"),
    "BtoKstarlow_AL0_err_noq2": ("B_Ks", "18_1_5"),
    "BtoKstarlow_AR0_err_noq2": ("B_Ks", "18_1_6"),
    "BtoKstarlow_At_err_noq2": ("B_Ks", "18_1_7"),
    "BtoKstarlow_AS_err_noq2": ("B_Ks", "18_1_8"),

    # --- B -> K* low-q2 q2-dependent piece ---
    "BtoKstarlow_ALperp_err_q2": ("B_Ks", "18_2_1"),
    "BtoKstarlow_ARperp_err_q2": ("B_Ks", "18_2_2"),
    "BtoKstarlow_ALpar_err_q2": ("B_Ks", "18_2_3"),
    "BtoKstarlow_ARpar_err_q2": ("B_Ks", "18_2_4"),
    "BtoKstarlow_AL0_err_q2": ("B_Ks", "18_2_5"),
    "BtoKstarlow_AR0_err_q2": ("B_Ks", "18_2_6"),
    "BtoKstarlow_At_err_q2": ("B_Ks", "18_2_7"),
    "BtoKstarlow_AS_err_q2": ("B_Ks", "18_2_8"),

    # --- B -> K* high-q2 piece ---
    "BtoKstarhigh_ALperp_err": ("B_Ks", "18_3_1"),
    "BtoKstarhigh_ARperp_err": ("B_Ks", "18_3_2"),
    "BtoKstarhigh_ALpar_err": ("B_Ks", "18_3_3"),
    "BtoKstarhigh_ARpar_err": ("B_Ks", "18_3_4"),
    "BtoKstarhigh_AL0_err": ("B_Ks", "18_3_5"),
    "BtoKstarhigh_AR0_err": ("B_Ks", "18_3_6"),
    "BtoKstarhigh_At_err": ("B_Ks", "18_3_7"),
    "BtoKstarhigh_AS_err": ("B_Ks", "18_3_8"),

    # --- B -> K low-q2 non-factorizable nuisances ---
    "BtoKlow_FV_err_noq2": ("B_K", "18_1_1"),
    "BtoKlow_FA_err_noq2": ("B_K", "18_1_2"),
    "BtoKlow_FS_err_noq2": ("B_K", "18_1_3"),
    "BtoKlow_FP_err_noq2": ("B_K", "18_1_4"),

    # --- B -> K low-q2 q2-dependent piece ---
    "BtoKlow_FV_err_q2": ("B_K", "18_2_1"),
    "BtoKlow_FA_err_q2": ("B_K", "18_2_2"),
    "BtoKlow_FS_err_q2": ("B_K", "18_2_3"),
    "BtoKlow_FP_err_q2": ("B_K", "18_2_4"),

    # --- B -> K high-q2 piece ---
    "BtoKhigh_FV_err": ("B_K", "18_3_1"),
    "BtoKhigh_FA_err": ("B_K", "18_3_2"),
    "BtoKhigh_FS_err": ("B_K", "18_3_3"),
    "BtoKhigh_FP_err": ("B_K", "18_3_4"),

     # --- Bs -> phi low-q2 non-factorizable nuisances ---
    "Bstophilow_ALperp_err_noq2": ("B_phi", "18_1_1"),
    "Bstophilow_ARperp_err_noq2": ("B_phi", "18_1_2"),
    "Bstophilow_ALpar_err_noq2": ("B_phi", "18_1_3"),
    "Bstophilow_ARpar_err_noq2": ("B_phi", "18_1_4"),
    "Bstophilow_AL0_err_noq2": ("B_phi", "18_1_5"),
    "Bstophilow_AR0_err_noq2": ("B_phi", "18_1_6"),
    "Bstophilow_At_err_noq2": ("B_phi", "18_1_7"),
    "Bstophilow_AS_err_noq2": ("B_phi", "18_1_8"),

    # --- Bs -> phi low-q2 q2-dependent piece ---
    "Bstophilow_ALperp_err_q2": ("B_phi", "18_2_1"),
    "Bstophilow_ARperp_err_q2": ("B_phi", "18_2_2"),
    "Bstophilow_ALpar_err_q2": ("B_phi", "18_2_3"),
    "Bstophilow_ARpar_err_q2": ("B_phi", "18_2_4"),
    "Bstophilow_AL0_err_q2": ("B_phi", "18_2_5"),
    "Bstophilow_AR0_err_q2": ("B_phi", "18_2_6"),
    "Bstophilow_At_err_q2": ("B_phi", "18_2_7"),
    "Bstophilow_AS_err_q2": ("B_phi", "18_2_8"),

    # --- Bs -> phi high-q2 piece ---
    "Bstophihigh_ALperp_err": ("B_phi", "18_3_1"),
    "Bstophihigh_ARperp_err": ("B_phi", "18_3_2"),
    "Bstophihigh_ALpar_err": ("B_phi", "18_3_3"),
    "Bstophihigh_ARpar_err": ("B_phi", "18_3_4"),
    "Bstophihigh_AL0_err": ("B_phi", "18_3_5"),
    "Bstophihigh_AR0_err": ("B_phi", "18_3_6"),
    "Bstophihigh_At_err": ("B_phi", "18_3_7"),
    "Bstophihigh_AS_err": ("B_phi", "18_3_8"),
}


# ---------------------------------------------------------------------------
# 2) One legacy nuisance -> several new entries
# ---------------------------------------------------------------------------

LEGACY_TO_NEW_MULTI_MAP: Dict[str, NewTargets] = {
    # same f_B used for neutral and charged B in the new file
    "f_B": [("FCONST", "511_1"), ("FCONST", "521_1")],

    # same f_K used for K0 and K+
    "f_K": [("FCONST", "311_1"), ("FCONST", "321_1")],

    # lambda_B-like input duplicated across channels
    "lambda_Bp": [("B_K", "13"), ("B_Ks", "13")],

    # a common lambda_Bs-type input for Bs -> phi
    "lambda_Bsp": [("B_phi", "13")],

    # HO notation in the legacy file has one shared a0 for both gperp / gplus sectors
    # depending on the exact convention. Keep both as candidate targets.
    "a0_HO_gpp_LbLll": [("Lb_L", "1_3_0"), ("Lb_L", "1_6_0")],

    # same comment for htilde "pp" in the pasted legacy notation
    "a0_HO_htildepp_LbLll": [("Lb_L", "1_7_0"), ("Lb_L", "1_8_0")],
}


# ---------------------------------------------------------------------------
# 3) Old detailed nuisances compressed/shared in the new format
# ---------------------------------------------------------------------------

COMPRESSED_OR_SHARED_MAP: Dict[str, NewTarget] = {
    
}


# ---------------------------------------------------------------------------
# 5) Legacy IDs seen in the old file but not mapped with enough confidence
# ---------------------------------------------------------------------------

UNMATCHED_OR_NOT_EXPOSED_IN_NEW: List[str] = [
    # Scales
    "log_mu_spec_lambda_h_mass_b",
    "log_muK_1GeV",
    "log_mu_W_mass_W",
    "log_mu_b_mass_b",
    
    # K_L,S > ll
    "Aterm_mu_KLmumu",
    "Iterm_mu_KSmumu",
    "chi_gg_Mrho",

    # K > pi l l (NYI in HI)
    "err_Pc_Xlambda_Kppipnunu",
    "deltaPcu_Kppipnunu",
    "KLpill_Ce_Cdir",
    "KLpill_Cmu_Cdir",
    "KLpill_Ce_Cint",
    "KLpill_Cmu_Cint",
    "KLpill_Ce_Cmix",
    "KLpill_Cmu_Cmix",
    "KLpill_Cmu_CPC",
    "KLpill_abs_aS",
    "KLpill_Ce_Sgg",
    "KLpill_Cmu_Sgg",
]


def get_new_nuisance_target(
    legacy_id: str,
) -> Union[NewTarget, NewTargets, None]:
    """
    Resolve a legacy nuisance id to the new-format target.

    Resolution order:
      1) exact one-to-one map
      2) exact one-to-many map
      3) compressed/shared map
      4) approximate map
      5) None if not mapped
    """
    if legacy_id in LEGACY_TO_NEW_NUISANCE_MAP:
        return LEGACY_TO_NEW_NUISANCE_MAP[legacy_id]
    if legacy_id in LEGACY_TO_NEW_MULTI_MAP:
        return LEGACY_TO_NEW_MULTI_MAP[legacy_id]
    if legacy_id in COMPRESSED_OR_SHARED_MAP:
        return COMPRESSED_OR_SHARED_MAP[legacy_id]
    return None


def is_unmatched_legacy_nuisance(legacy_id: str) -> bool:
    return legacy_id in UNMATCHED_OR_NOT_EXPOSED_IN_NEW


ALL_MAPPED_KEYS = set(LEGACY_TO_NEW_NUISANCE_MAP)
ALL_MAPPED_KEYS |= set(LEGACY_TO_NEW_MULTI_MAP)
ALL_MAPPED_KEYS |= set(COMPRESSED_OR_SHARED_MAP)
