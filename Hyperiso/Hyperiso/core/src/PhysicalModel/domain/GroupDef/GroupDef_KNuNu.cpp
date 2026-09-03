#include "GroupDefinition.h"
#include "KNuNuWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
const GroupDefinition& KNuNu() {
    static const GroupDefinition def = [] {
        GroupDefinition d;
        d.id = GroupMapper::to_id(WGroup::KNuNu);
        d.members = WCoefMapper::get_group(WGroup::KNuNu);

        std::map<QCDOrder, CGS> orders;
        CGS lo;
        lo.sources = {{ParameterType::WILSON, {
            MATCHING_BLOCK_PLACEHOLDER,
            "WPARAM_MATCH_SM",
            "WPARAM_RUN_SM"
        }}};
        lo.func = &KNuNuCoefficientGroup::base_1_LO_calculation;
        orders[QCDOrder::LO] = lo;

        CGS nlo = lo;
        nlo.func = &KNuNuCoefficientGroup::base_1_NLO_calculation;
        orders[QCDOrder::NLO] = nlo;

        CGS nnlo = lo;
        nnlo.func = &KNuNuCoefficientGroup::base_1_NNLO_calculation;
        orders[QCDOrder::NNLO] = nnlo;

        // Semileptonic colour-singlet vector operators: identity evolution at
        // the level of the full coefficient in the present WET basis.
        d.sources.emplace(WilsonBasis::B_STANDARD, std::move(orders));
        return d;
    }();
    return def;
}
}
