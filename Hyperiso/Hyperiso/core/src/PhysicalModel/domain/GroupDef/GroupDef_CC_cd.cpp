#include "GroupDefinition.h"
#include "DChargedCurrentsWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& CC_cd() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::CC_cd);
            d.members = WCoefMapper::get_group(WGroup::CC_cd);

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER } }
            };
            lo.func = &DdlnuCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
