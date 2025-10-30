#include "GroupDefinition.h"
#include "ChargedCurrentsWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& CC_bc() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::CC_bc);
            d.members = { WCoef::C_V1_bc, WCoef::C_V2_bc, WCoef::C_S1_bc, WCoef::C_S2_bc, WCoef::C_T_bc };

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER } }
            };
            lo.func = &BclnuCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
