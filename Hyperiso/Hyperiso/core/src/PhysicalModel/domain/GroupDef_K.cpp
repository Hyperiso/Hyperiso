#include "GroupDefinition.h"
#include "KWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& K() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::K);
            d.members = { WCoef::CK9, WCoef::CK10, WCoef::CKQ1, WCoef::CKQ2, WCoef::CPK9,
                          WCoef::CPK10, WCoef::CPKQ1, WCoef::CPKQ2, WCoef::CK_L};
            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER} }
            };
            lo.func = &KCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
