// main.cpp
#include <cassert>
#include <iostream>
#include <vector>

#include "General.h"          // LhaID
#include "GeneralEnum.h"      // enums (WCoef, Observables, Decays, etc.)
#include "registry_init.hpp"  // init_all_builtins()
#include "mapper_hub.hpp"     // unified access helpers
#include "decay_graph.hpp"    // DecayGraph::instance()

int main() {
    std::cout << "== Dynamic registry demo ==\n";

    // 1) Initialize all builtin mappers (one-time, idempotent).
    init_all_builtins();

    // -----------------------------------------------------------------------------
    // 2) WCOEF: string -> Id -> canonical, external roundtrip
    // -----------------------------------------------------------------------------
    {
        // Case-insensitive lookup (aliases supported by the runtime registry)
        auto id = WCoefMapper::id_of("c7");
        std::cout << "[WCOEF] id_of(\"c7\") = " << WCoefMapper::str(id) << "\n";
        assert(WCoefMapper::str(id) == "C7");

        // External key roundtrip (FLHA-style pair<int,int>)
        auto ext = WCoefMapper::external_of(id);
        assert(ext.has_value());
        auto back = WCoefMapper::from_external(*ext);
        assert(back && WCoefMapper::str(*back) == "C7");
        std::cout << "[WCOEF] FLHA roundtrip ok\n";
    }

    // Register a custom WCoef with an external key
    {
        bool ok = WCoefMapper::register_custom("CNEW", {"c_new"}, std::make_pair(999, 888));
        assert(ok);

        auto id  = WCoefMapper::id_of("CNEW");
        auto ext = WCoefMapper::external_of(id);

        assert(ext && ext->first == 999 && ext->second == 888);
        std::cout << "[WCOEF] custom CNEW registered with FLHA (999,888)\n";
    }

    // -----------------------------------------------------------------------------
    // 3) OBSERVABLES: enum -> Id -> LhaID -> back (+ legacy enum_elt(LhaID))
    // -----------------------------------------------------------------------------
    {
        // Start from a builtin enum
        auto e  = Observables::BR_BS_MUMU;
        auto id = ObservableMapper::to_id(e);

        // Get the external FLHA-like key
        auto lha = ObservableMapper::flha_of(id);
        assert(lha.has_value());

        // External back to runtime Id
        auto id2 = ObservableMapper::from_flha(*lha);
        assert(id2 && ObservableMapper::str(*id2) == ObservableMapper::str(id));

        // Legacy: directly to builtin enum from external
        auto e2 = ObservableMapper::enum_elt(*lha);
        assert(e2 == e);

        std::cout << "[OBS] LhaID roundtrip ok for " << ObservableMapper::str(id) << "\n";
    }

    // Register a custom Observable with LhaID and attach it to a decay
    {
        LhaID custom_key(777, 1, 2, 13, -13); // arbitrary example

        // Attach to builtin decay B__l_l; this also creates the runtime link in DecayGraph
        bool ok = ObservableMapper::register_custom(
            "MY_OBS",                 // canonical name
            {"myObs"},                // aliases
            custom_key,               // external key (LhaID)
            Decays::B__l_l           // parent decay (enum overload)
        );
        assert(ok);

        // Lookup by alias
        auto id = ObservableMapper::id_of("myObs");
        assert(ObservableMapper::str(id) == "MY_OBS");

        // External must match
        auto k = ObservableMapper::flha_of(id);
        assert(k && k->to_string() == custom_key.to_string());

        // DecayGraph now contains the link (only our custom one is required here)
        auto d  = DecayMapper::to_id(Decays::B__l_l);
        auto os = DecayGraph::instance().observables_of(d);

        bool found = false;
        for (const auto& x : os) {
            if (x.str() == "MY_OBS") { found = true; break; }
        }
        assert(found);
        std::cout << "[OBS] custom MY_OBS linked to decay B__l_l\n";
    }

    // -----------------------------------------------------------------------------
    // 4) DECAY: builtin usage (external key optional; none by default)
    // -----------------------------------------------------------------------------
    {
        auto d = DecayMapper::to_id(Decays::B__D_l_nu);
        std::cout << "[DECAY] canonical = " << DecayMapper::str(d) << "\n";

        // If you decide later to bind an external code to a builtin decay:
        // DecayMapper::set_external(d, LhaID(521));         // example (arbitrary)
        // auto ext = DecayMapper::external_of(d);           // would return that LhaID
        // assert(ext && static_cast<long>(*ext) == 521);    // depending on your LhaID implementation
    }

    // -----------------------------------------------------------------------------
    // 5) GROUP: block naming helpers (MATCHING / HADRONIC + basis)
    // -----------------------------------------------------------------------------
    {
        auto m = GroupMapper::str(WGroup::BCC, ScaleType::MATCHING);
        auto h = GroupMapper::str(WGroup::BCC, ScaleType::HADRONIC, WilsonBasis::B_TRADITIONAL);
        std::cout << "[GROUP] MATCHING=" << m << "  HADRONIC(TRAD)=" << h << "\n";
    }

    // -----------------------------------------------------------------------------
    // 6) QCD ORDER: string → runtime id → builtin enum
    // -----------------------------------------------------------------------------
    {
        // Generic API: id_of("NLO") gives a runtime id; enum_of(id) recovers the builtin enum
        auto qid = OrderMapper::id_of("NLO");
        auto q   = OrderMapper::enum_of(qid);
        assert(q && *q == QCDOrder::NLO);
        std::cout << "[ORDER] string->id->enum OK (NLO)\n";
    }

    // -----------------------------------------------------------------------------
    // 7) Unified access via Hub (MapperKind) for simple tasks
    // -----------------------------------------------------------------------------
    {
        // Canonical name from user input without knowing which concrete mapper at call site
        std::cout << "[HUB] canonical(WCoef, \"c7\") = "
                  << canonical_of(MapperKind::WCoef, "c7") << "\n";

        // Register a custom WCoef through the hub helper with an external key
        bool ok = register_custom_wcoef("CNEW2", {"cnew2"}, {1000, 1001});
        assert(ok);
        std::cout << "[HUB] registered custom WCoef CNEW2 via hub helper\n";

        // List a few items from a mapper via Hub
        auto obs_list = list_all(MapperKind::Observable);
        std::cout << "[HUB] first observable = " << (obs_list.empty() ? "<empty>" : obs_list.front()) << "\n";
    }

    std::cout << "✅ Demo OK\n";
    return 0;
}
