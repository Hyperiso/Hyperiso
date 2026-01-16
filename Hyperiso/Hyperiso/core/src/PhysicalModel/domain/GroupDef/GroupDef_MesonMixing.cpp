#include "GroupDefinition.h"
#include "MesonMixingWilsonGroup.h"

using CGS = CoefficientGroupSources;

static void Setup_Mixing_RunningMatrices([[maybe_unused]] const BuildContext& ctx, [[maybe_unused]] CoefficientGroup& grp) {
    // TODO : compose ETA_POWS_MIXING, UM_MATRIX_5, UM_MATRIX_4 avec ctx.adapters.iblock_c->compose_block(...)
}

namespace GroupDefinitions {
    const GroupDefinition& MesonMixing() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::MESON_MIXING);

            d.members = WCoefMapper::get_group(WGroup::MESON_MIXING);

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
