#include "MesonMixingWilsonGroup.h"

using MMRP = MesonMixingRunningParameters;

MesonMixingCoefficientGroup::MesonMixingCoefficientGroup(WilsonGroupAdapterConfig adapters) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::MESON_MIXING);

}

std::shared_ptr<CoefficientGroup> MesonMixingCoefficientGroup::clone() const
{
    return std::make_shared<MesonMixingCoefficientGroup>(*this);
}


//Fake !!
void MesonMixingCoefficientGroup::init_running_parameter_blocks() {
    // WilsonParamComposer composer;

    LOG_DEBUG("Init running matrices blocks of Meson Mixing Coefficient group");
	std::unordered_map<ParameterType, std::vector<std::string>> eta_powers_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM"}}};

    auto eta_powers_func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta5 = src.get_val("WPARAM_RUN_SM",2);
        double eta4 = src.get_val("WPARAM_RUN_SM",3);

		for (int i = 0; i < MMRP::n_pows; ++i) {
            dep_block->store_or_assign(LhaID(1, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS_MIXING", LhaID(1, i)}, std::pow(eta5, (MMRP::ai)[i]), 0., 0.));
            dep_block->store_or_assign(LhaID(2, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS_MIXING", LhaID(2, i)}, std::pow(eta4, (MMRP::bi)[i]), 0., 0.));
		}
    };

    std::unordered_map<ParameterType, std::vector<std::string>> mtx_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM", "ETA_POWS_MIXING"}}};

    auto U_5_func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta_5 = src.get_val("WPARAM_RUN_SM",2);
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "UM_MATRIX_5", LhaID(n, k, l)};
        };

        using MMRP = MesonMixingRunningParameters;

        double U0, U1;

        // V
        double eta_V = src.get_val("ETA_POWS_MIXING",LhaID(1, 0));
        U0 = MMRP::a0_V_5 * eta_V;
        U1 = (MMRP::a1_V_5 + MMRP::b_V_5 * eta_5) * eta_V;
        dep_block->store_or_assign(LhaID(0, 0, 0), std::make_shared<Parameter>(pid(0, 0, 0), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(0, 5, 5), std::make_shared<Parameter>(pid(0, 5, 5), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 0, 0), std::make_shared<Parameter>(pid(1, 0, 0), U1, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 5, 5), std::make_shared<Parameter>(pid(1, 5, 5), U1, 0., 0.));

        // LR
        std::array<double, 2> eta_LR {
            src.get_val("ETA_POWS_MIXING",LhaID(1, 1)),
            src.get_val("ETA_POWS_MIXING",LhaID(1, 2)),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) {
                    U0 += MMRP::a0_LR_5[i][j][k] * eta_LR[k];    
                    U1 += (MMRP::a1_LR_5[i][j][k] + MMRP::b_LR_5[i][j][k] * eta_5) * eta_LR[k];
                }
                dep_block->store_or_assign(LhaID(0, 1+i, 1+j), std::make_shared<Parameter>(pid(0, 1+i, 1+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 1+i, 1+j), std::make_shared<Parameter>(pid(1, 1+i, 1+j), U1, 0., 0.));
            }
        }

        // S
		std::array<double, 2> eta_S {
            src.get_val("ETA_POWS_MIXING",LhaID(1, 3)),
            src.get_val("ETA_POWS_MIXING",LhaID(1, 4)),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) {
                    U0 += MMRP::a0_S_5[i][j][k] * eta_S[k];    
                    U1 += (MMRP::a1_S_5[i][j][k] + MMRP::b_S_5[i][j][k] * eta_5) * eta_S[k];
                }
                dep_block->store_or_assign(LhaID(0, 3+i, 3+j), std::make_shared<Parameter>(pid(0, 3+i, 3+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 3+i, 3+j), std::make_shared<Parameter>(pid(1, 3+i, 3+j), U1, 0., 0.));
                dep_block->store_or_assign(LhaID(0, 6+i, 6+j), std::make_shared<Parameter>(pid(0, 6+i, 6+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 6+i, 6+j), std::make_shared<Parameter>(pid(1, 6+i, 6+j), U1, 0., 0.));
            }
        }
    };

    auto U_4_func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta_5 = src.get_val("WPARAM_RUN_SM",2);
        double eta_4 = src.get_val("WPARAM_RUN_SM",3);
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "UM_MATRIX_4", LhaID(n, k, l)};
        };

        using MMRP = MesonMixingRunningParameters;

        double U0, U1;

        // V
        double eta_5_V = src.get_val("ETA_POWS_MIXING",LhaID(1, 0));
        double eta_4_V = src.get_val("ETA_POWS_MIXING",LhaID(2, 0));
        U0 = MMRP::a0_V_4 * eta_4_V * eta_5_V;
        U1 = (MMRP::a1_V_4 + MMRP::b_V_4 * eta_4 + MMRP::c_V_4 * eta_4 * eta_5) * eta_4_V * eta_5_V;
        dep_block->store_or_assign(LhaID(0, 0, 0), std::make_shared<Parameter>(pid(0, 0, 0), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(0, 5, 5), std::make_shared<Parameter>(pid(0, 5, 5), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 0, 0), std::make_shared<Parameter>(pid(1, 0, 0), U1, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 5, 5), std::make_shared<Parameter>(pid(1, 5, 5), U1, 0., 0.));

        // LR
        std::array<double, 2> eta_5_LR {
            src.get_val("ETA_POWS_MIXING",LhaID(1, 1)),
            src.get_val("ETA_POWS_MIXING",LhaID(1, 2)),
        };

        std::array<double, 2> eta_4_LR {
            src.get_val("ETA_POWS_MIXING",LhaID(2, 1)),
            src.get_val("ETA_POWS_MIXING",LhaID(2, 2)),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) for (int l = 0; l < 2; ++l) {
                    U0 += MMRP::a0_LR_4[i][j][k][l] * eta_4_LR[k] * eta_5_LR[l];    
                    U1 += (MMRP::a1_LR_4[i][j][k][l] + MMRP::b_LR_4[i][j][k][l] * eta_4 + MMRP::c_LR_4[i][j][k][l] * eta_4 * eta_5) * eta_4_LR[k] * eta_5_LR[l];
                }
                dep_block->store_or_assign(LhaID(0, 1+i, 1+j), std::make_shared<Parameter>(pid(0, 1+i, 1+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 1+i, 1+j), std::make_shared<Parameter>(pid(1, 1+i, 1+j), U1, 0., 0.));
            }
        }

        // S
		std::array<double, 2> eta_5_S {
            src.get_val("ETA_POWS_MIXING",LhaID(1, 3)),
            src.get_val("ETA_POWS_MIXING",LhaID(1, 4)),
        };

        std::array<double, 2> eta_4_S {
            src.get_val("ETA_POWS_MIXING",LhaID(2, 3)),
            src.get_val("ETA_POWS_MIXING",LhaID(2, 4)),
        };

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                U0 = U1 = 0;
                for (int k = 0; k < 2; ++k) for (int l = 0; l < 2; ++l) {
                    U0 += MMRP::a0_S_4[i][j][k][l] * eta_4_S[k] * eta_5_S[l];    
                    U1 += (MMRP::a1_S_4[i][j][k][l] + MMRP::b_S_4[i][j][k][l] * eta_4 + MMRP::c_S_4[i][j][k][l] * eta_4 * eta_5) * eta_4_S[k] * eta_5_S[l];
                }
                dep_block->store_or_assign(LhaID(0, 3+i, 3+j), std::make_shared<Parameter>(pid(0, 3+i, 3+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 3+i, 3+j), std::make_shared<Parameter>(pid(1, 3+i, 3+j), U1, 0., 0.));
                dep_block->store_or_assign(LhaID(0, 6+i, 6+j), std::make_shared<Parameter>(pid(0, 6+i, 6+j), U0, 0., 0.));
                dep_block->store_or_assign(LhaID(1, 6+i, 6+j), std::make_shared<Parameter>(pid(1, 6+i, 6+j), U1, 0., 0.));
            }
        }
    };

    adapters.iblock_c->compose_block("ETA_POWS_MIXING", eta_powers_src, eta_powers_func);
    adapters.iblock_c->compose_block("UM_MATRIX_5", mtx_src, U_5_func);
    adapters.iblock_c->compose_block("UM_MATRIX_4", mtx_src, U_4_func);
}

std::unordered_map<WCoefId, scalar_t> 
MesonMixingCoefficientGroup::base_1_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
    const BlockSrc& src)
{
    int n_f_final = QCDHelper::get_nf(src.get_val("B_SCALE",1));
    std::string src_block = n_f_final < 5 ? "UM_MATRIX_4" : "UM_MATRIX_5"; 

    std::array<complex_t, 32> Ci_match_BMU = {};
    std::array<complex_t, 8> Ci_match_temp = {};
    auto ids = WCoefMapper::get_group(WGroup::MESON_MIXING);
    for (size_t n = 0; n < 4; n++) {
        for (size_t k = 0; k < 8; k++) {
            Ci_match_temp[k] = coef_matching.at(QCDOrder::LO).at(WCoefMapper::to_id(ids[8 * n + k]));
        }
        Ci_match_temp = MMRP::change_basis(Ci_match_temp, MMRP::SUSY_to_BMU);
        for (size_t k = 0; k < 8; k++) {
            Ci_match_BMU[8 * n + k] = Ci_match_temp[k];
        }
        Ci_match_temp = MMRP::change_basis(Ci_match_temp, MMRP::BMU_to_SUSY);
    }    

    double fact = 4. * PI / src.get_val("WPARAM_RUN_SM",1);

    std::array<complex_t, 32> Ci_run {};

    for (size_t k = 0; k < 8; k++) {
        for (size_t l = 0; l < 8; l++) {
            double U0 = src.raw().at(src_block)->contains(LhaID(0, k, l)) ? src.raw().at(src_block)->retrieve(LhaID(0, k, l))->get_val() : scalar_t();
            double U1 = src.raw().at(src_block)->contains(LhaID(1, k, l)) ? src.raw().at(src_block)->retrieve(LhaID(1, k, l))->get_val() : scalar_t();
            double U = U0 + fact * U1;
            Ci_run[k] += U * Ci_match_BMU[l];
            Ci_run[k + 8] += U * Ci_match_BMU[l + 8];
            Ci_run[k + 16] += U * Ci_match_BMU[l + 16];
            Ci_run[k + 24] += U * Ci_match_BMU[l + 24];
        }
    }

    std::unordered_map<WCoefId, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 32; k++) {
        Ci_run_map[WCoefMapper::to_id(ids[k])] = Ci_run[k];
    }

    return Ci_run_map;
}