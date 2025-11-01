#include "GroupDefinition.h"
#include "BWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& BPrime() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::BPrime);
            d.members = { WCoef::CP1, WCoef::CP2, WCoef::CP3, WCoef::CP4, WCoef::CP5,
                          WCoef::CP6, WCoef::CP7, WCoef::CP8, WCoef::CP9, WCoef::CP10,
                          WCoef::CPQ1, WCoef::CPQ2 };
            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
            };
            lo.func = &BPrimeCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            return d;
        }();
        return def;
    }
}
