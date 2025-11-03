#include <cassert>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "Include.h"
#include "Wilson.h"
#include "WilsonGroup.h"

#include "GroupDefinition.h"
#include "CoefficientGroupBuilder.h"
#include "WilsonCoefficientRegistry.h"

// ====== Dummies / Spies ======
class DummyBoolAPI    : public ICoreAPI<bool>      { public: bool get() override { return false; } };
class DummyStringAPI  : public ICoreAPI<std::string>{ public: std::string get() override { return "SM"; } };
class DummyPathAPI    : public ICoreAPI<fs::path>  { public: fs::path get() override { return fs::path("/dev/null"); } };

class SpyProxy : public IParameterProxy<std::string, LhaID> {
public:
    mutable std::string last_block;
    mutable LhaID       last_id{0};
    double ret = 5.0;
    scalar_t operator()(const std::string& b, const LhaID& id) const override {
        last_block=b; last_id=id; return ret;
    }
    bool   exist(const std::string&, const LhaID&) const override { return true; }
    double get_scale(const std::string&) const override { return 0.0; }
};

struct ParamRecord { ParamId target; std::unordered_set<ParamId> sources; };

class SpyComposer : public IBlockComposer {
public:
    std::vector<ParamRecord> params;

    void compose_block(const std::string&,
                       const std::unordered_map<ParameterType, std::vector<std::string>>&,
                       const DepUpdateFunc&) override {}

    void compose_parameter(const ParamId& id,
                           const std::unordered_set<ParamId>& src,
                           const DepParamUpdateFunc&) override {
        params.push_back({id, src});
    }

    void remove_block(const std::string&) override {}
    void update(const std::string&) override {}
    void remove_all_composed_blocks() override {}
};

// ====== Registry : makers pour B (SM builtin + Marty si besoin) ======
extern void register_B(CoefficientRegistry& reg);

// -------- helpers “mini-manager” (hadronic only) --------

// dé-sérialise un LhaID wilson -> (WCoef, (QCDOrder, ContributionType))
static std::pair<WCoef, std::pair<QCDOrder, ContributionType>> lha_wilson_deserialize(LhaID id) {
    auto parts = id.get_parts(); // [w1, w2, order, contrib]
    auto w_id = std::make_pair<int, int>(parts[0], parts[1]);
    if (!WCoefMapper::inverse_flha_mapping().contains(w_id)) {
        LOG_ERROR("ValueError", "bad lha id for wilson conversion");
    }
    WCoef coef = WCoefMapper::inverse_flha_mapping().at(w_id);
    QCDOrder order = parts[2] ? ((parts[2] -1) ? QCDOrder::NNLO : QCDOrder::NLO) : QCDOrder::LO;
    ContributionType part = parts[3] ? parts[3] -1 ? ContributionType::TOTAL : ContributionType::BSM : ContributionType::SM;
    return {coef, {order, part}};
}

// merge des sources LO..ord pour un basis donné
static void fill_sources_for_group(std::shared_ptr<CoefficientGroup> grp,
                                   QCDOrder ord,
                                   std::unordered_map<ParameterType, std::vector<std::string>>& src,
                                   WilsonBasis basis)
{
    auto merge = [&](QCDOrder o){
        for (const auto& [ptype, vec] : grp->get_sources(o, basis)) {
            auto& dst = src[ptype];
            dst.insert(dst.end(), vec.begin(), vec.end());
        }
    };
    if (ord >= QCDOrder::LO)  merge(QCDOrder::LO);
    if (ord >= QCDOrder::NLO) merge(QCDOrder::NLO);
    if (ord >= QCDOrder::NNLO)merge(QCDOrder::NNLO);
}

// composition du bloc hadronic pour un basis (copie light du Manager)
static void init_group_hadronic_no_manager(std::shared_ptr<CoefficientGroup> grp,
                                           WilsonGroupAdapterConfig& adapters,
                                           QCDOrder ord,
                                           WilsonBasis basis)
{
    std::unordered_map<ParameterType, std::vector<std::string>> src;
    fill_sources_for_group(grp, ord, src, basis);

    // fonctions LO/NLO/NNLO (si absentes, la map contiendra des std::function vides -> à gérer)
    std::map<QCDOrder, std::function<std::unordered_map<WCoef, scalar_t>(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>&,
        const BlockSrc&
    )>> funcs = {
        {QCDOrder::LO,   grp->get_func(QCDOrder::LO,   basis)},
        {QCDOrder::NLO,  grp->get_func(QCDOrder::NLO,  basis)},
        {QCDOrder::NNLO, grp->get_func(QCDOrder::NNLO, basis)}
    };

    std::string matching_block_name = grp->get_matching_storage_block();
    std::string hadronic_block_name = GroupMapper::str(grp->get_group_id(), ScaleType::HADRONIC, basis);

    auto func = [matching_block_name, ord, funcs, grp, basis]
        (const BlockSrc& src,
         std::shared_ptr<DependentBlock> dep_block)
    {
        // 1) lire tous les matching coefs déposés dans le bloc MATCHING
        std::map<LhaID, std::shared_ptr<Parameter>> matching_coeff = src.raw().at(matching_block_name)->getItems();
        std::unordered_map<ContributionType, std::unordered_map<QCDOrder, std::unordered_map<WCoef, scalar_t>>> matching_map;

        for (auto& kv : matching_coeff) {
            auto des = lha_wilson_deserialize(kv.first);
            const WCoef& wcoef = des.first;
            const QCDOrder& order = des.second.first;
            const ContributionType& contrib = des.second.second;
            matching_map[contrib][order][wcoef] = kv.second->get_val();
        }

        // 2) exécuter les fonctions LO..ord pour SM/BSM/TOTAL
        std::unordered_map<ContributionType, std::unordered_map<QCDOrder,std::unordered_map<WCoef, scalar_t>>> res;
        auto call_if = [&](QCDOrder o) {
            auto f = funcs.at(o);
            if (!f) return std::unordered_map<WCoef, scalar_t>{};
            return f(matching_map[ContributionType::SM], BlockSrc(src)); // on calcule par contrib juste après
        };

        for (auto contri : {ContributionType::SM, ContributionType::BSM, ContributionType::TOTAL}) {
            switch (ord) {
            case QCDOrder::NNLO:
                if (funcs.at(QCDOrder::NNLO)) res[contri][QCDOrder::NNLO] = funcs.at(QCDOrder::NNLO)(matching_map[contri], BlockSrc(src));
                [[fallthrough]];
            case QCDOrder::NLO:
                if (funcs.at(QCDOrder::NLO))  res[contri][QCDOrder::NLO]  = funcs.at(QCDOrder::NLO)(matching_map[contri], BlockSrc(src));
                [[fallthrough]];
            case QCDOrder::LO:
                if (funcs.at(QCDOrder::LO))   res[contri][QCDOrder::LO]   = funcs.at(QCDOrder::LO)(matching_map[contri], BlockSrc(src));
                break;
            default: break;
            }
        }

        // 3) write dans le bloc hadronic
        for (auto& [c_type, order_map] : res) {
            for (auto& [order, coef_map] : order_map) {
                for (auto& [coef_id, coef_val] : coef_map) {
                    LhaID coef_lha = WCoefMapper::flha_full(coef_id, order, c_type);
                    ParamId pid {
                        ParameterType::WILSON,
                        GroupMapper::str(grp->get_group_id(), ScaleType::HADRONIC, basis),
                        coef_lha
                    };
                    dep_block->store_or_assign(pid.code, std::make_shared<Parameter>(pid, coef_val, 0., (int)c_type));
                }
            }
        }
    };

    adapters.iblock_c->compose_block(hadronic_block_name, src, func);
}

int main() {
    std::cout << "== INTEGRATION (no manager) ==\n";

    // --- Adapters ---
    auto proxy   = std::make_shared<SpyProxy>();
    auto comp    = std::make_shared<SpyComposer>();
    auto use_mty = std::make_shared<DummyBoolAPI>(); // false ⇒ Backend::Builtin
    auto mname   = std::make_shared<DummyStringAPI>();
    auto mpath   = std::make_shared<DummyPathAPI>();
    WilsonGroupAdapterConfig adapters(proxy, comp, use_mty, mname, mpath);

    // --- Registry + Builder ---
    CoefficientRegistry reg;
    register_B(reg); // enregistre C1..C10 pour SM Builtin (+ Marty SM si tu l’as mis)
    CoefficientGroupBuilder builder(reg);

    // --- Construit le groupe B via GroupDefinitions + Registry (SANS init) ---
    BuildContext ctx{
        .adapters = adapters,
        .model    = Model::SM,
        .backend  = Backend::Builtin,
        .contrib  = ContributionType::SM,
        .group_id = GroupMapper::to_id(WGroup::B)
    };
    auto grp = builder.build(ctx);
    assert(grp->size() == 10);

    // --- Matching seulement (comme avant, sans manager)
    grp->init(QCDOrder::LO);

    // vérif dépendance C7@LO sur WPARAM_MATCH_SM
    const auto match_blk = GroupMapper::str(WGroup::B, ScaleType::MATCHING); // même nom qu’avant
    auto c7 = grp->at("C7");
    auto c7_lo = c7->get_lhaid(QCDOrder::LO);

    bool saw_c7_lo_dep = false;
    for (const auto& r : comp->params) {
        if (r.target == ParamId{match_blk, c7_lo}) {
            ParamId need{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2,1)};
            if (r.sources.count(need) == 1) { saw_c7_lo_dep = true; break; }
        }
    }
    assert(saw_c7_lo_dep);

    // --- “mini-manager” pour hadronic (sinon, pas de bloc hadronic composé)
    init_group_hadronic_no_manager(grp, adapters, QCDOrder::LO, WilsonBasis::B_STANDARD);

    // --- Lecture run
    proxy->ret = 9.0;
    auto vr = grp->get_running_coefficient("C10", "LO", ContributionType::SM, WilsonBasis::B_STANDARD);
    assert(std::abs(vr.real() - 9.0) < 1e-12);

    auto had_blk = GroupMapper::str(WGroup::B, ScaleType::HADRONIC, WilsonBasis::B_STANDARD);
    assert(proxy->last_block == had_blk);
    assert(proxy->last_id == grp->at("C10")->id(QCDOrder::LO, ContributionType::SM));

    std::cout << "✅ INTEGRATION OK (no manager)\n";
    return 0;
}
