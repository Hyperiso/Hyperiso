#include "GroupDefinition.h"
#include "BWilsonGroup.h"

using CGS = CoefficientGroupSources;

// Hook SUSY (portage de BScalarCoefficientGroup_susy::set_base_1_LO)
static std::unordered_map<WCoef, scalar_t>
BScalar_SUSY_Base1_LO_calculation(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>& coef_matching,
    const std::unordered_map<std::string, std::shared_ptr<Block>>& src)
{
    // Membres du groupe (ex: CQ1, CQ2)
    const auto ids = WCoefMapper::get_group(WGroup::BScalar);

    // Accès aux matching LO (robuste)
    const auto itLO = coef_matching.find(QCDOrder::LO);
    const auto* matchLO = (itLO != coef_matching.end()) ? &itLO->second : nullptr;

    auto getM = [&](WCoef c) -> scalar_t {
        if (!matchLO) return scalar_t(0);
        auto it = matchLO->find(c);
        return (it != matchLO->end()) ? it->second : scalar_t(0);
    };

    // Paramètres de running (identiques à la "normale")
    const double eta    = src.at("WPARAM_RUN_SM")->retrieve(2)->get_val();
    const double beta_0 = src.at("WPARAM_SI_SM")->retrieve(5)->get_val(); // TODO: passer aux vrais params QCD si dispo

    const double fact = std::pow(eta, -4.0 / beta_0);

    // Running LO : même logique que la version normale
    std::unordered_map<WCoef, scalar_t> out;
    out.reserve(ids.size());
    for (auto c : ids) {
        out[c] = fact * getM(c);
    }

    // ====== Point d’extension pour corrections SUSY (facultatif, à réintégrer plus tard) ======
    // Exemple de squelette si tu veux ajouter un delta sur CQ1/CQ2 :
    // {
    //     // Lis ici GAUGE/HMIX/... via src.at("...")->retrieve(...)
    //     // calcule delta_CQ1, delta_CQ2 (complex_t/scalar_t)
    //     // puis :
    //     // out[WCoef::CQ1] += delta_CQ1;
    //     // out[WCoef::CQ2] += delta_CQ2;
    // }

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
    lo.func = &BScalar_SUSY_Base1_LO_calculation;
    m[QCDOrder::LO] = lo;
    
    CGS nlo;
    nlo.sources = {
        // calque la définition par défaut : pas besoin de MATCH_SM ici
        { ParameterType::WILSON, { GroupMapper::str(ctx.group_id, ScaleType::MATCHING),
                                   "WPARAM_RUN_SM", "WPARAM_SI_SM" } }
    };
    nlo.func = &BScalarCoefficientGroup::base_1_NLO_calculation;
    m[QCDOrder::NLO] = nlo;

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
