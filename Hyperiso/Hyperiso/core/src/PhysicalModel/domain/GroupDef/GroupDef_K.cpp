#include "GroupDefinition.h"
#include "KWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& K() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::K);
            d.members = WCoefMapper::get_group(WGroup::K);
            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER} }
            };
            lo.func = &KCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;


            CGS nlo = lo;
            nlo.func = &KCoefficientGroup::base_1_NLO_calculation;

            CGS nnlo = lo;
            nnlo.func = &KCoefficientGroup::base_1_NNLO_calculation;

            m[QCDOrder::NLO]  = nlo;
            m[QCDOrder::NNLO] = nnlo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
