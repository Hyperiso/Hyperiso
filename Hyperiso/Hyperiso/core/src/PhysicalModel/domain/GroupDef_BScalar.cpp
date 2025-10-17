#include "GroupDefinition.h"
#include "BWilsonGroup.h"

using CGS = CoefficientGroupSources;

// Hook SUSY (portage de BScalarCoefficientGroup_susy::set_base_1_LO)
static std::unordered_map<WCoef, scalar_t>
BScalar_SUSY_LO_running(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src)
{
    auto getM = [&](WCoef c)->scalar_t {
        auto itO = coef_matching.find(QCDOrder::LO);
        if (itO == coef_matching.end()) return scalar_t(0);
        auto itC = itO->second.find(c);
        return itC == itO->second.end() ? scalar_t(0) : itC->second;
    };

    // paramètres nécessaires (identiques à ton code)
    double eta    = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    double beta_0 = src.at("WPARAM_SI_SM")->retrieve(5)->get_val();

    auto powfac = std::pow(eta, -4.0 / beta_0);

    std::unordered_map<WCoef, scalar_t> out;
    out[WCoef::CQ1] = getM(WCoef::CQ1) * powfac;
    out[WCoef::CQ2] = getM(WCoef::CQ2) * powfac;

    // NOTE : tu peux réintégrer ici, pas à pas, le reste de ta correction SUSY
    // en t'appuyant sur 'src' pour GAUGE/HMIX/... et sur getM(...) pour les matching coefs.
    return out;
}

// Hook SUSY : on *override* LO pour BScalar (basis standard) + on ajoute les sources BSM
static void Setup_BScalar_SUSY_Base1_LO(const BuildContext& ctx, CoefficientGroup& grp) {
    // re-déclare *spécifiquement* la LO pour SUSY avec les dépendances adéquates
    std::map<QCDOrder, CGS> m;

    CGS lo;
    lo.sources = {
        // IMPORTANT : matching block sans basis (le Manager l'écrit ainsi)
        { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
                                   "WPARAM_RUN_SM", "WPARAM_SI_SM", "WPARAM_MATCH_SM" } },
        { ParameterType::SM,     { "MASS", "SMINPUTS" } },
        { ParameterType::BSM,    { "GAUGE", "HMIX", "STOPMIX", "UMIX", "VMIX", "NMAMIX", "NMHMIX", "MASS" } }
    };
    lo.func = &BScalar_SUSY_LO_running;
    m[QCDOrder::LO] = lo;

    // On laisse NLO tel quel (hérité de la def par défaut) si tu veux le garder :
    // sinon, tu peux également surcharger NLO ici de la même manière.
    grp.add_sources(WilsonBasis::B_STANDARD, std::move(m));


}

namespace GroupDefinitions {
    const GroupDefinition& BScalar() {
        static const GroupDefinition def = []{
            GroupDefinition d;
            d.id = GroupMapper::to_id(WGroup::BScalar);
            d.members = { WCoef::CQ1, WCoef::CQ2 };

            std::map<QCDOrder, CGS> m;
            CGS lo;
            lo.sources = {
                { ParameterType::WILSON, { MATCHING_BLOCK_PLACEHOLDER, "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
            };
            lo.func = &BScalarCoefficientGroup::base_1_LO_calculation;
            m[QCDOrder::LO] = lo;

            CGS nlo = lo; nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
            m[QCDOrder::NLO] = nlo;

            d.sources.emplace(WilsonBasis::B_STANDARD, std::move(m));

            d.setup[Model::SUSY].push_back(&Setup_BScalar_SUSY_Base1_LO);
            return d;
        }();
        return def;
    }
}
