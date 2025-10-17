#include "GroupDefinition.h"
#include "MesonMixingWilsonGroup.h"

using CGS = CoefficientGroupSources;

static void Setup_Mixing_RunningMatrices(const BuildContext& ctx, CoefficientGroup& grp) {
    (void)grp;
    // === recopier ici le contenu de ton init_running_parameter_blocks() ===
    // compose ETA_POWS_MIXING, UM_MATRIX_5, UM_MATRIX_4 avec ctx.adapters.iblock_c->compose_block(...)
    // (je te laisse coller tes lambdas U_5_func / U_4_func telles quelles)
}

namespace GroupDefinitions {
    const GroupDefinition& MesonMixing() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::MESON_MIXING);

            d.members = {
                WCoef::C_BD_1,  WCoef::CT_BD_1, WCoef::C_BD_2,  WCoef::CT_BD_2,  WCoef::C_BD_3,  WCoef::CT_BD_3,  WCoef::C_BD_4,  WCoef::C_BD_5,
                WCoef::C_BS_1,  WCoef::CT_BS_1, WCoef::C_BS_2,  WCoef::CT_BS_2,  WCoef::C_BS_3,  WCoef::CT_BS_3,  WCoef::C_BS_4,  WCoef::C_BS_5,
                WCoef::C_SD_1,  WCoef::CT_SD_1, WCoef::C_SD_2,  WCoef::CT_SD_2,  WCoef::C_SD_3,  WCoef::CT_SD_3,  WCoef::C_SD_4,  WCoef::C_SD_5,
                WCoef::C_CU_1,  WCoef::CT_CU_1, WCoef::C_CU_2,  WCoef::CT_CU_2,  WCoef::C_CU_3,  WCoef::CT_CU_3,  WCoef::C_CU_4,  WCoef::C_CU_5
            };

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "UM_MATRIX_5", "UM_MATRIX_4", "B_SCALE" } }
            };
            lo.func = &MesonMixingCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));

            d.setup[Model::SM].push_back(&Setup_Mixing_RunningMatrices);
            d.setup[Model::SUSY].push_back(&Setup_Mixing_RunningMatrices);
            d.setup[Model::THDM].push_back(&Setup_Mixing_RunningMatrices);
            return d;
        }();
        return def;
    }
}
