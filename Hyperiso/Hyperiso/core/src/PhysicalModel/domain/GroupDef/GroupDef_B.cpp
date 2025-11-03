#include "GroupDefinition.h"
#include "BWilsonGroup.h"

using CGS = CoefficientGroupSources;

namespace GroupDefinitions {
    const GroupDefinition& B() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::B);
            d.members = WCoefMapper::get_group(WGroup::B);

            {
                std::map<QCDOrder, CGS> m;
                CGS lo;
                lo.sources = {
                    { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "WPARAM_MATCH_SM", "B_SCALE", "U_MATRIX" } },
                    { ParameterType::SM,     { "SMINPUTS", "MASS" } }
                };
                lo.func = &BCoefficientGroup::base_1_LO_calculation;
                m[QCDOrder::LO] = lo;

                CGS nlo = lo;  nlo.func  = &BCoefficientGroup::base_1_NLO_calculation;
                m[QCDOrder::NLO] = nlo;

                CGS nnlo = lo; nnlo.func = &BCoefficientGroup::base_1_NNLO_calculation;
                m[QCDOrder::NNLO] = nnlo;

                d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));
            }
            {
                std::map<QCDOrder, CGS> m;
                CGS lo;
                lo.sources = {
                    { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "WPARAM_MATCH_SM", "V_MATRIX" } }
                };
                lo.func = &BCoefficientGroup::base_2_LO_calculation;
                m[QCDOrder::LO] = lo;

                CGS nlo = lo; nlo.func = &BCoefficientGroup::base_2_NLO_calculation;
                m[QCDOrder::NLO] = nlo;

                d.sources.emplace(WilsonBasis::B_TRADITIONAL, std::move(m));
            }
            return d;
        }();
        return def;
    }
}
