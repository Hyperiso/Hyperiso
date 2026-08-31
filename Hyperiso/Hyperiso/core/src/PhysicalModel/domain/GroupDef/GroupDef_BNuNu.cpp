#include "GroupDefinition.h"
#include "BNuNuWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
const GroupDefinition& BNuNu() {
    static const GroupDefinition def = [] {
        GroupDefinition d;
        d.id = GroupMapper::to_id(WGroup::BNuNu);
        d.members = WCoefMapper::get_group(WGroup::BNuNu);

        std::map<QCDOrder, CGS> orders;
        CGS lo;
        lo.sources = {{ParameterType::WILSON, {MATCHING_BLOCK_PLACEHOLDER}}};
        lo.func = &BNuNuCoefficientGroup::base_1_LO_calculation;
        orders[QCDOrder::LO] = lo;

        CGS nlo = lo;
        nlo.func = &BNuNuCoefficientGroup::base_1_NLO_calculation;
        orders[QCDOrder::NLO] = nlo;

        CGS nnlo = lo;
        nnlo.func = &BNuNuCoefficientGroup::base_1_NNLO_calculation;
        orders[QCDOrder::NNLO] = nnlo;

        // O_L/R are semileptonic colour-singlet vector operators.  No
        // non-trivial QCD mixing is applied in this WET basis.
        d.sources.emplace(WilsonBasis::B_STANDARD, std::move(orders));
        return d;
    }();
    return def;
}
}
