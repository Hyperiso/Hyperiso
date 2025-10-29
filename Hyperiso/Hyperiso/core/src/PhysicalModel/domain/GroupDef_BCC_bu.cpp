#include "GroupDefinition.h"
#include "ChargedCurrentsWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& BCC_bu() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::BCC_bu);
            d.members = { WCoef::C_V1_bu, WCoef::C_V2_bu, WCoef::C_S1_bu, WCoef::C_S2_bu, WCoef::C_T_bu };

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER } }
            };
            lo.func = &BulnuCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
