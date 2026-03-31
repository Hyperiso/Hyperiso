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
    "log_mu_W_mass_W": ("EW_SCALE", "1"),
    "log_mu_b_mass_b": ("B_SCALE", "1"),
    "log_muK_1GeV": ("K_SCALE", "1"),

    # --- inclusive B -> Xs ---
    "BR_BXclnu_exp": ("B_Xs", "2"),
    "mu_G2_bsg": ("B_Xs", "3"),
    "rho_D3_bsg": ("B_Xs", "4"),
    "rho_LS3_bsg": ("B_Xs", "5"),

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

    # --- Lambda_b -> Lambda ---
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
    # --- B -> K* low-q2 non-factorizable nuisances ---
    "BtoKstarlow_ALperp_err_noq2": ("B_Ks", "18_1"),
    "BtoKstarlow_ARperp_err_noq2": ("B_Ks", "18_1"),
    "BtoKstarlow_ALpar_err_noq2": ("B_Ks", "18_1"),
    "BtoKstarlow_ARpar_err_noq2": ("B_Ks", "18_1"),
    "BtoKstarlow_AL0_err_noq2": ("B_Ks", "18_1"),
    "BtoKstarlow_AR0_err_noq2": ("B_Ks", "18_1"),
    "BtoKstarlow_At_err_noq2": ("B_Ks", "18_1"),
    "BtoKstarlow_AS_err_noq2": ("B_Ks", "18_1"),

    # --- B -> K* low-q2 q2-dependent piece ---
    "BtoKstarlow_ALperp_err_q2": ("B_Ks", "18_2"),
    "BtoKstarlow_ARperp_err_q2": ("B_Ks", "18_2"),
    "BtoKstarlow_ALpar_err_q2": ("B_Ks", "18_2"),
    "BtoKstarlow_ARpar_err_q2": ("B_Ks", "18_2"),
    "BtoKstarlow_AL0_err_q2": ("B_Ks", "18_2"),
    "BtoKstarlow_AR0_err_q2": ("B_Ks", "18_2"),
    "BtoKstarlow_At_err_q2": ("B_Ks", "18_2"),
    "BtoKstarlow_AS_err_q2": ("B_Ks", "18_2"),

    # --- B -> K* high-q2 piece ---
    "BtoKstarhigh_ALperp_err": ("B_Ks", "18_3"),
    "BtoKstarhigh_ARperp_err": ("B_Ks", "18_3"),
    "BtoKstarhigh_ALpar_err": ("B_Ks", "18_3"),
    "BtoKstarhigh_ARpar_err": ("B_Ks", "18_3"),
    "BtoKstarhigh_AL0_err": ("B_Ks", "18_3"),
    "BtoKstarhigh_AR0_err": ("B_Ks", "18_3"),
    "BtoKstarhigh_At_err": ("B_Ks", "18_3"),
    "BtoKstarhigh_AS_err": ("B_Ks", "18_3"),
}


# ---------------------------------------------------------------------------
# 4) Approximate matches:
#    same role / same sector, but changed packaging or noticeably changed values
# ---------------------------------------------------------------------------

APPROXIMATE_MAP: Dict[str, NewTarget] = {
    # process-specific charm scale nuisance; closest target in the new inclusive block
    "mu_c_bsg": ("B_Xs", "7"),

    # generic scale nuisance reused across several channel blocks in the new layout
    "log_mu_spec_lambda_h_mass_b": ("B_Ks", "14"),

    # the new kaon blocks look like a compressed rewrite of the older KL->pi ll parameterization
    # these are placeholders only if you want a best-effort anchor:
    "deltaPcu_Kppipnunu": ("K_lnu", "1"),
    "log_muK_1GeV": ("K_ll", "4"),
}


# ---------------------------------------------------------------------------
# 5) Legacy IDs seen in the old file but not mapped with enough confidence
# ---------------------------------------------------------------------------

UNMATCHED_OR_NOT_EXPOSED_IN_NEW: List[str] = [
    # inclusive B -> X_s(ll) random nuisance placeholders
    "bsgamma_rand",
    "BRBXsmumu_lowq2_rand",
    "BRBXsmumu_highq2_rand",
    "BRBXsmumu_full_rand",
    "BRBXsee_lowq2_rand",
    "BRBXsee_highq2_rand",
    "BRBXsee_full_rand",

    # old charm-loop polynomial / alpha-expansion nuisances
    "real_alpha_perp0",
    "real_alpha_perp1",
    "real_alpha_perp2",
    "real_alpha_par0",
    "real_alpha_par1",
    "real_alpha_par2",
    "real_alpha_zero0",
    "real_alpha_zero1",
    "imag_alpha_perp0",
    "imag_alpha_perp1",
    "imag_alpha_perp2",
    "imag_alpha_par0",
    "imag_alpha_par1",
    "imag_alpha_par2",
    "imag_alpha_zero0",
    "imag_alpha_zero1",
    "DeltaC9_M1_q2bar",
    "r1_M1",
    "r2_M1",
    "DeltaC9_M2_q2bar",
    "r1_M2",
    "r2_M2",
    "DeltaC9_M3_q2bar",
    "r1_M3",
    "r2_M3",

    # B -> K low/high non-factorizable nuisance basis not explicitly visible in the new excerpt
    "BtoKlow_FV_err_noq2",
    "BtoKlow_FA_err_noq2",
    "BtoKlow_FS_err_noq2",
    "BtoKlow_FP_err_noq2",
    "BtoKlow_FV_err_q2",
    "BtoKlow_FA_err_q2",
    "BtoKlow_FS_err_q2",
    "BtoKlow_FP_err_q2",
    "BtoKhigh_FV_err",
    "BtoKhigh_FA_err",
    "BtoKhigh_FS_err",
    "BtoKhigh_FP_err",

    # Bs -> phi low/high detailed nuisance basis not explicitly exposed in the new excerpt
    "Bstophilow_ALperp_err_noq2",
    "Bstophilow_ARperp_err_noq2",
    "Bstophilow_ALpar_err_noq2",
    "Bstophilow_ARpar_err_noq2",
    "Bstophilow_AL0_err_noq2",
    "Bstophilow_AR0_err_noq2",
    "Bstophilow_At_err_noq2",
    "Bstophilow_AS_err_noq2",
    "Bstophilow_ALperp_err_q2",
    "Bstophilow_ARperp_err_q2",
    "Bstophilow_ALpar_err_q2",
    "Bstophilow_ARpar_err_q2",
    "Bstophilow_AL0_err_q2",
    "Bstophilow_AR0_err_q2",
    "Bstophilow_At_err_q2",
    "Bstophilow_AS_err_q2",
    "Bstophihigh_ALperp_err",
    "Bstophihigh_ARperp_err",
    "Bstophihigh_ALpar_err",
    "Bstophihigh_ARpar_err",
    "Bstophihigh_AL0_err",
    "Bstophihigh_AR0_err",
    "Bstophihigh_At_err",
    "Bstophihigh_AS_err",

    # lifetime of Lambda_b not exposed in the pasted new blocks
    "life_Lb",

    # old kaon-ll and KL->pi ll coefficients: appear repackaged in new K_ll / K_pi,
    # but not one-to-one from the pasted excerpt
    "err_Pc_Xlambda_Kppipnunu",
    "Aterm_mu_KLmumu",
    "chi_gg_Mrho",
    "Iterm_mu_KSmumu",
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

    # alternative fit/model families from the legacy file; the pasted new JSON only shows one active family
    "a0A0_BKstar_GRvDV_BSZ",
    "a1A0_BKstar_GRvDV_BSZ",
    "a2A0_BKstar_GRvDV_BSZ",
    "a0A1_BKstar_GRvDV_BSZ",
    "a1A1_BKstar_GRvDV_BSZ",
    "a2A1_BKstar_GRvDV_BSZ",
    "a1A12_BKstar_GRvDV_BSZ",
    "a2A12_BKstar_GRvDV_BSZ",
    "a0V_BKstar_GRvDV_BSZ",
    "a1V_BKstar_GRvDV_BSZ",
    "a2V_BKstar_GRvDV_BSZ",
    "a0T1_BKstar_GRvDV_BSZ",
    "a1T1_BKstar_GRvDV_BSZ",
    "a2T1_BKstar_GRvDV_BSZ",
    "a1T2_BKstar_GRvDV_BSZ",
    "a2T2_BKstar_GRvDV_BSZ",
    "a0T23_BKstar_GRvDV_BSZ",
    "a1T23_BKstar_GRvDV_BSZ",
    "a2T23_BKstar_GRvDV_BSZ",

    "a0A0_BKstar_GKvD_LCSR_Lattice",
    "a1A0_BKstar_GKvD_LCSR_Lattice",
    "a2A0_BKstar_GKvD_LCSR_Lattice",
    "a0A1_BKstar_GKvD_LCSR_Lattice",
    "a1A1_BKstar_GKvD_LCSR_Lattice",
    "a2A1_BKstar_GKvD_LCSR_Lattice",
    "a1A12_BKstar_GKvD_LCSR_Lattice",
    "a2A12_BKstar_GKvD_LCSR_Lattice",
    "a0V_BKstar_GKvD_LCSR_Lattice",
    "a1V_BKstar_GKvD_LCSR_Lattice",
    "a2V_BKstar_GKvD_LCSR_Lattice",
    "a0T1_BKstar_GKvD_LCSR_Lattice",
    "a1T1_BKstar_GKvD_LCSR_Lattice",
    "a2T1_BKstar_GKvD_LCSR_Lattice",
    "a1T2_BKstar_GKvD_LCSR_Lattice",
    "a2T2_BKstar_GKvD_LCSR_Lattice",
    "a0T23_BKstar_GKvD_LCSR_Lattice",
    "a1T23_BKstar_GKvD_LCSR_Lattice",
    "a2T23_BKstar_GKvD_LCSR_Lattice",

    "a0A0_BKstar_GKvD_LCSRonly",
    "a1A0_BKstar_GKvD_LCSRonly",
    "a2A0_BKstar_GKvD_LCSRonly",
    "a0A1_BKstar_GKvD_LCSRonly",
    "a1A1_BKstar_GKvD_LCSRonly",
    "a2A1_BKstar_GKvD_LCSRonly",
    "a1A12_BKstar_GKvD_LCSRonly",
    "a2A12_BKstar_GKvD_LCSRonly",
    "a0V_BKstar_GKvD_LCSRonly",
    "a1V_BKstar_GKvD_LCSRonly",
    "a2V_BKstar_GKvD_LCSRonly",
    "a0T1_BKstar_GKvD_LCSRonly",
    "a1T1_BKstar_GKvD_LCSRonly",
    "a2T1_BKstar_GKvD_LCSRonly",
    "a1T2_BKstar_GKvD_LCSRonly",
    "a2T2_BKstar_GKvD_LCSRonly",
    "a0T23_BKstar_GKvD_LCSRonly",
    "a1T23_BKstar_GKvD_LCSRonly",
    "a2T23_BKstar_GKvD_LCSRonly",

    "a0A0_Bsphi_GRvDV_BSZ",
    "a1A0_Bsphi_GRvDV_BSZ",
    "a2A0_Bsphi_GRvDV_BSZ",
    "a0A1_Bsphi_GRvDV_BSZ",
    "a1A1_Bsphi_GRvDV_BSZ",
    "a2A1_Bsphi_GRvDV_BSZ",
    "a1A12_Bsphi_GRvDV_BSZ",
    "a2A12_Bsphi_GRvDV_BSZ",
    "a0V_Bsphi_GRvDV_BSZ",
    "a1V_Bsphi_GRvDV_BSZ",
    "a2V_Bsphi_GRvDV_BSZ",
    "a0T1_Bsphi_GRvDV_BSZ",
    "a1T1_Bsphi_GRvDV_BSZ",
    "a2T1_Bsphi_GRvDV_BSZ",
    "a1T2_Bsphi_GRvDV_BSZ",
    "a2T2_Bsphi_GRvDV_BSZ",
    "a0T23_Bsphi_GRvDV_BSZ",
    "a1T23_Bsphi_GRvDV_BSZ",
    "a2T23_Bsphi_GRvDV_BSZ",

    "a0fp_BK_GRvDV_BSZ",
    "a1fp_BK_GRvDV_BSZ",
    "a2fp_BK_GRvDV_BSZ",
    "a1f0_BK_GRvDV_BSZ",
    "a2f0_BK_GRvDV_BSZ",
    "a0fT_BK_GRvDV_BSZ",
    "a1fT_BK_GRvDV_BSZ",
    "a2fT_BK_GRvDV_BSZ",

    "a1f0_BK_GKvD_LCSR_Lattice",
    "a2f0_BK_GKvD_LCSR_Lattice",
    "a0fT_BK_GKvD_LCSR_Lattice",
    "a1fT_BK_GKvD_LCSR_Lattice",
    "a2fT_BK_GKvD_LCSR_Lattice",
    "a0fp_BK_GKvD_LCSR_Lattice",
    "a1fp_BK_GKvD_LCSR_Lattice",
    "a2fp_BK_GKvD_LCSR_Lattice",

    "a1fp_BK_GKvD_LCSRonly",
    "a2fp_BK_GKvD_LCSRonly",
    "a0fp_BK_GKvD_LCSRonly",
    "a1f0_BK_GKvD_LCSRonly",
    "a2f0_BK_GKvD_LCSRonly",
    "a1fT_BK_GKvD_LCSRonly",
    "a2fT_BK_GKvD_LCSRonly",
    "a0fT_BK_GKvD_LCSRonly",

    "a0V_BKstar_BSZ_LCSRonly",
    "a1V_BKstar_BSZ_LCSRonly",
    "a2V_BKstar_BSZ_LCSRonly",
    "a0T1_BKstar_BSZ_LCSRonly",
    "a1T1_BKstar_BSZ_LCSRonly",
    "a2T1_BKstar_BSZ_LCSRonly",
    "a0T2_BKstar_BSZ_LCSRonly",
    "a1T2_BKstar_BSZ_LCSRonly",
    "a2T2_BKstar_BSZ_LCSRonly",
    "a0T23_BKstar_BSZ_LCSRonly",
    "a1T23_BKstar_BSZ_LCSRonly",
    "a2T23_BKstar_BSZ_LCSRonly",
    "a0A0_BKstar_BSZ_LCSRonly",
    "a1A0_BKstar_BSZ_LCSRonly",
    "a2A0_BKstar_BSZ_LCSRonly",
    "a0A1_BKstar_BSZ_LCSRonly",
    "a1A1_BKstar_BSZ_LCSRonly",
    "a2A1_BKstar_BSZ_LCSRonly",
    "a0A12_BKstar_BSZ_LCSRonly",
    "a1A12_BKstar_BSZ_LCSRonly",
    "a2A12_BKstar_BSZ_LCSRonly",

    "a0V_Bsphi_BSZ_LCSRonly",
    "a1V_Bsphi_BSZ_LCSRonly",
    "a2V_Bsphi_BSZ_LCSRonly",
    "a0T1_Bsphi_BSZ_LCSRonly",
    "a1T1_Bsphi_BSZ_LCSRonly",
    "a2T1_Bsphi_BSZ_LCSRonly",
    "a0T2_Bsphi_BSZ_LCSRonly",
    "a1T2_Bsphi_BSZ_LCSRonly",
    "a2T2_Bsphi_BSZ_LCSRonly",
    "a0T23_Bsphi_BSZ_LCSRonly",
    "a1T23_Bsphi_BSZ_LCSRonly",
    "a2T23_Bsphi_BSZ_LCSRonly",
    "a0A0_Bsphi_BSZ_LCSRonly",
    "a1A0_Bsphi_BSZ_LCSRonly",
    "a2A0_Bsphi_BSZ_LCSRonly",
    "a0A1_Bsphi_BSZ_LCSRonly",
    "a1A1_Bsphi_BSZ_LCSRonly",
    "a2A1_Bsphi_BSZ_LCSRonly",
    "a0A12_Bsphi_BSZ_LCSRonly",
    "a1A12_Bsphi_BSZ_LCSRonly",
    "a2A12_Bsphi_BSZ_LCSRonly",

    "a0V_BKstar_HLMW",
    "a1V_BKstar_HLMW",
    "a0A0_BKstar_HLMW",
    "a1A0_BKstar_HLMW",
    "a0A1_BKstar_HLMW",
    "a1A1_BKstar_HLMW",
    "a0A12_BKstar_HLMW",
    "a1A12_BKstar_HLMW",
    "a0T1_BKstar_HLMW",
    "a1T1_BKstar_HLMW",
    "a0T2_BKstar_HLMW",
    "a1T2_BKstar_HLMW",
    "a0T23_BKstar_HLMW",
    "a1T23_BKstar_HLMW",
    "systErr_BKstar_HLMW",

    "a0V_Bsphi_HLMW",
    "a1V_Bsphi_HLMW",
    "a0A0_Bsphi_HLMW",
    "a1A0_Bsphi_HLMW",
    "a0A1_Bsphi_HLMW",
    "a1A1_Bsphi_HLMW",
    "a0A12_Bsphi_HLMW",
    "a1A12_Bsphi_HLMW",
    "a0T1_Bsphi_HLMW",
    "a1T1_Bsphi_HLMW",
    "a0T2_Bsphi_HLMW",
    "a1T2_Bsphi_HLMW",
    "a0T23_Bsphi_HLMW",
    "a1T23_Bsphi_HLMW",
    "systErr_Bsphi_HLMW",

    "a0fp_BK_FLAG24",
    "a1fp_BK_FLAG24",
    "a2fp_BK_FLAG24",
    "a0f0_BK_FLAG24",
    "a1f0_BK_FLAG24",
    "a0fT_BK_FLAG24",
    "a1fT_BK_FLAG24",
    "a2fT_BK_FLAG24",

    "a0f0_BK_HPQCD22",
    "a1f0_BK_HPQCD22",
    "a2f0_BK_HPQCD22",
    "a1fp_BK_HPQCD22",
    "a2fp_BK_HPQCD22",
    "a0fT_BK_HPQCD22",
    "a1fT_BK_HPQCD22",
    "a2fT_BK_HPQCD22",

    "f0fp_BK_KR",
    "b1fp_BK_KR",
    "f0fT_BK_KR",
    "b1fT_BK_KR",
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
    if legacy_id in APPROXIMATE_MAP:
        return APPROXIMATE_MAP[legacy_id]
    return None


def is_unmatched_legacy_nuisance(legacy_id: str) -> bool:
    return legacy_id in UNMATCHED_OR_NOT_EXPOSED_IN_NEW


ALL_MAPPED_KEYS = set(LEGACY_TO_NEW_NUISANCE_MAP)
ALL_MAPPED_KEYS |= set(LEGACY_TO_NEW_MULTI_MAP)
ALL_MAPPED_KEYS |= set(COMPRESSED_OR_SHARED_MAP)
ALL_MAPPED_KEYS |= set(APPROXIMATE_MAP)
