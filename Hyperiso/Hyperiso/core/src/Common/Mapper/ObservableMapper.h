#ifndef OBSERVABLE_MAPPER_H
#define OBSERVABLE_MAPPER_H

#include "GeneralEnum.h"
#include "EnumMapper.h"
#include "General.h"

class ObservableMapper : public EnumMapperBase<Observables, ObservableMapper> {
public:
    static const std::map<Observables, std::string>& mapping() {
        static const std::map<Observables, std::string> m = {
            {Observables::BR_BS_MUMU, "BR_Bs__mu_mu"},
            {Observables::BR_BS_MUMU_UNTAG, "BRuntag_Bs__mu_mu"},
            {Observables::BR_BD_MUMU, "BR_Bd__mu_mu"},
            {Observables::BR_BU_TAU_NU, "BR_Bu__tau_nu"},
            {Observables::R_TAU_NU, "R_tau_nu"},
            {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, "IA_B__K*_gamma"},
            {Observables::BR_B_XS_GAMMA, "BR_B__Xs_gamma"},
            {Observables::BR_B__D_TAU_NU, "BR_B__D_tau_nu"},
            {Observables::A_FB_B__D_TAU_NU, "A_FB_B__D_tau_nu"},
            {Observables::P_TAU_B__D_TAU_NU, "P_tau_B__D_tau_nu"},
            {Observables::R_D, "R_D"},
            {Observables::BR_B__DSTAR_TAU_NU, "BR_B__D*_tau_nu"},
            {Observables::A_FB_B__DSTAR_TAU_NU, "A_FB_B__D*_tau_nu"},
            {Observables::P_TAU_B__DSTAR_TAU_NU, "P_tau_B__D*_tau_nu"},
            {Observables::P_D_B__DSTAR_TAU_NU, "P_D*_B__D*_tau_nu"},
            {Observables::R_DSTAR, "R_D*"},
            {Observables::BR_B__Xs_mu_mu__LOW_Q2, "BR_B__Xs_mu_mu__[1-6]"},
            {Observables::BR_B__Xs_mu_mu__HIGH_Q2, "BR_B__Xs_mu_mu__[>14.4]"},
            {Observables::BR_B__Xs_tau_tau__HIGH_Q2, "BR_B__Xs_tau_tau__[>14.4]"},
            {Observables::PHI_D, "phi_d"},
            {Observables::DELTA_M_BD, "Delta_M_Bd"},
            {Observables::PHI_S, "phi_s"},
            {Observables::DELTA_M_BS, "Delta_M_Bs"},
            {Observables::A_FS, "a_fs"},
            {Observables::DELTA_M_K, "Delta_M_K"},
            {Observables::ABS_EPSILON_K, "|epsilon_K|"},
            {Observables::X_D, "x_D"},
            {Observables::TEST_B__KS_L_L, "test_B_Ks_ll"},
        };
        return m;
    }

    static const std::map<std::string, Observables>& inverse_mapping() {
        static const std::map<std::string, Observables> inv = invert_map(mapping());
        return inv;
    }

    static const std::map<Observables, LhaID>& flha_mapping() {
        static const std::map<Observables, LhaID> m = {
            {Observables::BR_BS_MUMU, LhaID(531, 1, 2, 13, -13)},
            {Observables::BR_BS_MUMU_UNTAG, LhaID(531, 15, 2, 13, -13)},
            {Observables::BR_BD_MUMU, LhaID(511, 1, 2, 13, -13)},
            {Observables::BR_BU_TAU_NU, LhaID(521, 1, 2, -15, 16)},
            {Observables::R_TAU_NU, LhaID(521, 2, 2, -15, 16)},
            {Observables::ISOSPIN_ASYMMETRY_B_KSTAR_GAMMA, LhaID(521, 4, 2, 313, 22)},
            {Observables::BR_B_XS_GAMMA, LhaID(5, 1, 2, 3, 22)},
            {Observables::BR_B__D_TAU_NU, LhaID(521, 1, 3, 421, -15, 16)},
            {Observables::A_FB_B__D_TAU_NU, LhaID(521, 5, 3, 421, -15, 16)},
            {Observables::P_TAU_B__D_TAU_NU, LhaID(521, 82, 3, 421, -15, 16)},
            {Observables::R_D, LhaID(521, 11, 3, 421, -15, 16)},
            {Observables::BR_B__DSTAR_TAU_NU, LhaID(521, 1, 3, 423, -15, 16)},
            {Observables::A_FB_B__DSTAR_TAU_NU, LhaID(521, 5, 3, 423, -15, 16)},
            {Observables::P_TAU_B__DSTAR_TAU_NU, LhaID(521, 82, 3, 423, -15, 16)},
            {Observables::P_D_B__DSTAR_TAU_NU, LhaID(521, 81, 3, 423, -15, 16)},
            {Observables::R_DSTAR, LhaID(521, 11, 3, 423, -15, 16)}
        };
        return m;
    }

    static const std::map<LhaID, Observables>& inverse_flha_mapping() {
        static const std::map<LhaID, Observables> inv = []() {
            std::map<LhaID, Observables> temp;
            for (const auto& [k, v] : flha_mapping()) temp[v] = k;
            return temp;
        }();
        return inv;
    }

    static LhaID flha(Observables obs) {
        return flha_mapping().at(obs);
    }

    static Observables enum_elt(LhaID id) {
        return inverse_flha_mapping().at(id);
    }
};



inline std::ostream& operator<<(std::ostream& os, const Observables& oid) {
    os << ObservableMapper::str(oid);
    return os;
};

#endif