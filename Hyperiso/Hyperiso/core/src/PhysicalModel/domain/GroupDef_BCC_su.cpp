#include "GroupDefinition.h"
#include "KChargedCurrentsWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& CC_su() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::CC_su);
            d.members = { WCoef::C_V1_su, WCoef::C_V2_su, WCoef::C_S1_su, WCoef::C_S2_su, WCoef::C_T_su };

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER } }
            };
            lo.func = &KulnuCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
