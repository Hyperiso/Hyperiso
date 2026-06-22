#include "GroupDefinition.h"
#include "MesonMixingWilsonGroup.h"
#include "WilsonParametersHelper.h"

using CGS = CoefficientGroupSources;

static void Setup_Mixing_RunningMatrices(const BuildContext& ctx, [[maybe_unused]] CoefficientGroup& grp) {
    WilsonParameterHelper::compose_meson_mixing_running_blocks(ctx.adapters.iblock_c);
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

            d.common_setup.push_back(&Setup_Mixing_RunningMatrices);
            return d;
        }();
        return def;
    }
}
