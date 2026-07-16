#include "GroupDefinition.h"
#include "DChargedCurrentsWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& CC_cs() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::CC_cs);
            d.members = WCoefMapper::get_group(WGroup::CC_cs);

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER } }
            };
            lo.func = &DslnuCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
