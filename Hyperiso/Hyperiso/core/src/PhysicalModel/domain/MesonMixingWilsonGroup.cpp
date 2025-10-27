#include "MesonMixingWilsonGroup.h"

using MMRP = MesonMixingRunningParameters;

MesonMixingCoefficientGroup::MesonMixingCoefficientGroup(WilsonGroupAdapterConfig adapters, bool force_sm) : CoefficientGroup(adapters) {
    this->id = GroupMapper::to_id(WGroup::MESON_MIXING);

}

std::shared_ptr<CoefficientGroup> MesonMixingCoefficientGroup::clone() const
{
    return std::make_shared<MesonMixingCoefficientGroup>(*this);
}

// void MesonMixingCoefficientGroup::init_sources() {
//     init_running_parameter_blocks();
//     std::map<QCDOrder,CoefficientGroupSources> grp_src;

//     grp_src[QCDOrder::LO].sources = {
//         {ParameterType::WILSON, {this->get_matching_storage_block(), "WPARAM_RUN_SM", "UM_MATRIX_5", "UM_MATRIX_4", "B_SCALE"}},
//     };

//     grp_src[QCDOrder::LO].func = base_1_LO_calculation;

//     this->sources.insert({WilsonBasis::B_STANDARD, grp_src});
// }

//Fake !!
void MesonMixingCoefficientGroup::init_running_parameter_blocks() {
    // WilsonParamComposer composer;

    LOG_DEBUG("Init running matrices blocks of Meson Mixing Coefficient group");
	std::unordered_map<ParameterType, std::vector<std::string>> eta_powers_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM"}}};

    auto eta_powers_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta5 = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        double eta4 = src.at("WPARAM_RUN_SM")->retrieve(3)->get_val();

		for (int i = 0; i < MMRP::n_pows; ++i) {
            dep_block->store_or_assign(LhaID(1, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS_MIXING", LhaID(1, i)}, std::pow(eta5, (MMRP::ai)[i]), 0., 0.));
            dep_block->store_or_assign(LhaID(2, i), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "ETA_POWS_MIXING", LhaID(2, i)}, std::pow(eta4, (MMRP::bi)[i]), 0., 0.));
		}
    };

    std::unordered_map<ParameterType, std::vector<std::string>> mtx_src = {{ParameterType::WILSON, {"WPARAM_RUN_SM", "ETA_POWS_MIXING"}}};

    auto U_5_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta_5 = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "UM_MATRIX_5", LhaID(n, k, l)};
        };

        using MMRP = MesonMixingRunningParameters;

        double U0, U1;

        // V
        double eta_V = src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 0))->get_val();
        U0 = MMRP::a0_V_5 * eta_V;
        U1 = (MMRP::a1_V_5 + MMRP::b_V_5 * eta_5) * eta_V;
        dep_block->store_or_assign(LhaID(0, 0, 0), std::make_shared<Parameter>(pid(0, 0, 0), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(0, 5, 5), std::make_shared<Parameter>(pid(0, 5, 5), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 0, 0), std::make_shared<Parameter>(pid(1, 0, 0), U1, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 5, 5), std::make_shared<Parameter>(pid(1, 5, 5), U1, 0., 0.));

        // LR
        std::array<double, 2> eta_LR {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 1))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 2))->get_val(),
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
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 3))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 4))->get_val(),
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

    auto U_4_func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double eta_5 = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
        double eta_4 = src.at("WPARAM_RUN_SM")->retrieve(3)->get_val();
        auto pid = [] (int n, int k, int l) {
            return ParamId{ParameterType::WILSON, "UM_MATRIX_4", LhaID(n, k, l)};
        };

        using MMRP = MesonMixingRunningParameters;

        double U0, U1;

        // V
        double eta_5_V = src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 0))->get_val();
        double eta_4_V = src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 0))->get_val();
        U0 = MMRP::a0_V_4 * eta_4_V * eta_5_V;
        U1 = (MMRP::a1_V_4 + MMRP::b_V_4 * eta_4 + MMRP::c_V_4 * eta_4 * eta_5) * eta_4_V * eta_5_V;
        dep_block->store_or_assign(LhaID(0, 0, 0), std::make_shared<Parameter>(pid(0, 0, 0), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(0, 5, 5), std::make_shared<Parameter>(pid(0, 5, 5), U0, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 0, 0), std::make_shared<Parameter>(pid(1, 0, 0), U1, 0., 0.));
        dep_block->store_or_assign(LhaID(1, 5, 5), std::make_shared<Parameter>(pid(1, 5, 5), U1, 0., 0.));

        // LR
        std::array<double, 2> eta_5_LR {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 1))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 2))->get_val(),
        };

        std::array<double, 2> eta_4_LR {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 1))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 2))->get_val(),
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
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 3))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(1, 4))->get_val(),
        };

        std::array<double, 2> eta_4_S {
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 3))->get_val(),
            src.at("ETA_POWS_MIXING")->retrieve(LhaID(2, 4))->get_val(),
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

std::unordered_map<WCoef, scalar_t> 
MesonMixingCoefficientGroup::base_1_LO_calculation (
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src)
{
    int n_f_final = QCDHelper::get_nf(src.at("B_SCALE")->retrieve(1)->get_val());
    std::string src_block = n_f_final < 5 ? "UM_MATRIX_4" : "UM_MATRIX_5"; 

    std::array<complex_t, 32> Ci_match_BMU = {};
    std::array<complex_t, 8> Ci_match_temp = {};
    auto ids = WCoefMapper::get_group(WGroup::MESON_MIXING);
    for (size_t n = 0; n < 4; n++) {
        for (size_t k = 0; k < 8; k++) {
            Ci_match_temp[k] = coef_matching.at(QCDOrder::LO).at(ids[8 * n + k]);
        }
        Ci_match_temp = MMRP::change_basis(Ci_match_temp, MMRP::SUSY_to_BMU);
        for (size_t k = 0; k < 8; k++) {
            Ci_match_BMU[8 * n + k] = Ci_match_temp[k];
        }
        Ci_match_temp = MMRP::change_basis(Ci_match_temp, MMRP::BMU_to_SUSY);
    }    

    double fact = 4 * PI / src.at("WPARAM_RUN_SM")->retrieve(1)->get_val();

    std::array<complex_t, 32> Ci_run {};

    for (size_t k = 0; k < 8; k++) {
        for (size_t l = 0; l < 8; l++) {
            double U0 = src.at(src_block)->contains(LhaID(0, k, l)) ? src.at(src_block)->retrieve(LhaID(0, k, l))->get_val() : scalar_t();
            double U1 = src.at(src_block)->contains(LhaID(1, k, l)) ? src.at(src_block)->retrieve(LhaID(1, k, l))->get_val() : scalar_t();
            double U = U0 + fact * U1;
            Ci_run[k] += U * Ci_match_BMU[l];
            Ci_run[k + 8] += U * Ci_match_BMU[l + 8];
            Ci_run[k + 16] += U * Ci_match_BMU[l + 16];
            Ci_run[k + 24] += U * Ci_match_BMU[l + 24];
        }
    }

    std::unordered_map<WCoef, scalar_t> Ci_run_map {};
    for (size_t k = 0; k < 32; k++) {
        Ci_run_map[ids[k]] = Ci_run[k];
    }

    return Ci_run_map;
}

// void MesonMixingCoefficientGroup::add_wilson_coefficients(bool force_sm) {
//     if (adapters.use_marty->get()) {
//         this->wilson_type = force_sm ? ContributionType::SM : ContributionType::TOTAL;
//         for (auto&& coeff : {"C_BD_1", "CT_BD_1", "C_BD_2", "CT_BD_2", "C_BD_3", "CT_BD_3", "C_BD_4", "C_BD_5", "C_BS_1", "CT_BS_1", "C_BS_2", "CT_BS_2", "C_BS_3", "CT_BS_3", "C_BS_4", "C_BS_5", "C_SD_1", "CT_SD_1", "C_SD_2", "CT_SD_2", "C_SD_3", "CT_SD_3", "C_SD_4", "C_SD_5", "C_CU_1", "CT_CU_1", "C_CU_2", "CT_CU_2", "C_CU_3", "CT_CU_3", "C_CU_4", "C_CU_5"}) {
//             std::string _name = force_sm ? "SM" : adapters.marty_model_name->get();
//             // std::string _name = force_sm ? "SM" : MartyModelNameAPI().get();
//             fs::path _path = force_sm ? adapters.sm_path : adapters.marty_model_path->get();
//             // fs::path _path = force_sm ? fs::path(std::string(project_assets_root.data())+"input_files/marty_model/sm.h") : MartyModelPathAPI().get();
//             std::string _block = GroupMapper::str(this->id, ScaleType::MATCHING);
//             LhaID _id = WCoefMapper::flha_full(WCoefMapper::enum_elt(coeff), QCDOrder::LO, this->get_type());
//             // this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(_id, _block, _name, _path)));
//             MartyWilsonConfig config {_name, _id, _block, _path, adapters.marty_proxy};
//             this->insert(std::make_pair(coeff, std::make_shared<MartyWilson>(config)));
//         }
//         return;
//     }

//     this->insert(std::make_pair("C_BD_1", std::make_shared<C_mix_bd_1>())); 
//     this->insert(std::make_pair("CT_BD_1", std::make_shared<C_mix_bd_1_tilde>())); 
//     this->insert(std::make_pair("C_BD_2", std::make_shared<C_mix_bd_2>()));
//     this->insert(std::make_pair("CT_BD_2", std::make_shared<C_mix_bd_2_tilde>()));  
//     this->insert(std::make_pair("C_BD_3", std::make_shared<C_mix_bd_3>())); 
//     this->insert(std::make_pair("CT_BD_3", std::make_shared<C_mix_bd_3_tilde>())); 
//     this->insert(std::make_pair("C_BD_4", std::make_shared<C_mix_bd_4>()));  
//     this->insert(std::make_pair("C_BD_5", std::make_shared<C_mix_bd_5>()));  

//     this->insert(std::make_pair("C_BS_1", std::make_shared<C_mix_bs_1>())); 
//     this->insert(std::make_pair("CT_BS_1", std::make_shared<C_mix_bs_1_tilde>()));
//     this->insert(std::make_pair("C_BS_2", std::make_shared<C_mix_bs_2>())); 
//     this->insert(std::make_pair("CT_BS_2", std::make_shared<C_mix_bs_2_tilde>())); 
//     this->insert(std::make_pair("C_BS_3", std::make_shared<C_mix_bs_3>())); 
//     this->insert(std::make_pair("CT_BS_3", std::make_shared<C_mix_bs_3_tilde>())); 
//     this->insert(std::make_pair("C_BS_4", std::make_shared<C_mix_bs_4>())); 
//     this->insert(std::make_pair("C_BS_5", std::make_shared<C_mix_bs_5>())); 

//     this->insert(std::make_pair("C_SD_1", std::make_shared<C_mix_sd_1>())); 
//     this->insert(std::make_pair("CT_SD_1", std::make_shared<C_mix_sd_1_tilde>())); 
//     this->insert(std::make_pair("C_SD_2", std::make_shared<C_mix_sd_2>())); 
//     this->insert(std::make_pair("CT_SD_2", std::make_shared<C_mix_sd_2_tilde>())); 
//     this->insert(std::make_pair("C_SD_3", std::make_shared<C_mix_sd_3>())); 
//     this->insert(std::make_pair("CT_SD_3", std::make_shared<C_mix_sd_3_tilde>())); 
//     this->insert(std::make_pair("C_SD_4", std::make_shared<C_mix_sd_4>()));
//     this->insert(std::make_pair("C_SD_5", std::make_shared<C_mix_sd_5>()));

//     this->insert(std::make_pair("C_CU_1", std::make_shared<C_mix_cu_1>()));
//     this->insert(std::make_pair("CT_CU_1", std::make_shared<C_mix_cu_1_tilde>()));
//     this->insert(std::make_pair("C_CU_2", std::make_shared<C_mix_cu_2>()));
//     this->insert(std::make_pair("CT_CU_2", std::make_shared<C_mix_cu_2_tilde>()));
//     this->insert(std::make_pair("C_CU_3", std::make_shared<C_mix_cu_3>()));
//     this->insert(std::make_pair("CT_CU_3", std::make_shared<C_mix_cu_3_tilde>()));
//     this->insert(std::make_pair("C_CU_4", std::make_shared<C_mix_cu_4>()));
//     this->insert(std::make_pair("C_CU_5", std::make_shared<C_mix_cu_5>()));
// }