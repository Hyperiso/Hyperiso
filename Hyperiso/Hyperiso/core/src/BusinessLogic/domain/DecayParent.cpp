#include "DecayParent.h"

// std::shared_ptr<WilsonInterface> DecayParent::get_wilsons(bool force_update) {
//     if (!wilson || force_update || true) {
//         wilson = std::make_shared<WilsonInterface>();
//         wilson->build(winfo.wgroups, winfo.matching_scale, winfo.hadronic_scale, winfo.order);
        
//         if (winfo.basis == BWilsonBasis::TRADITIONAL 
//             && std::find(winfo.wgroups.begin(), winfo.wgroups.end(), WGroup::B) != winfo.wgroups.end()) {
//             wilson->switchbasis(WGroup::B);
//         }
//     } else {
//         // for (auto &p : manager->getGroups()) {
//             // if (p.second->cacheChanged()) {
//                 // manager->update(p.first);
//             // }
//         // }
//     }

//     return wilson;
// }

//TODO : everything here, just for make it works
void DecayParent::set_order(QCDOrder new_order) {
    if (ObsUseMarty().get() && new_order > QCDOrder::LO) {
        LOG_WARN("Using MARTY defaults all calculations to LO in QCD.");
        new_order = QCDOrder::LO;
    }
    return;
    // if (winfo.order == new_order) {
    //     return;
    // }
    
    // if (winfo.order == QCDOrder::NONE) {
    //     if (new_order <= max_order && new_order != QCDOrder::NONE) {
    //         winfo.order = new_order;
    //     } else {
    //         winfo.order = max_order;
    //         LOG_WARN("Cannot set decay order to", OrderMapper::str(new_order));
    //     }
    // } else {
    //     LOG_WARN("QCD order for this decay has already been set.");
    // }
}

scalar_t DecayParent::compute_observable(Observables obs) {
    return roots.at(obs)->calculate();
}

size_t DecayParent::get_n_evals(Observables obs) {
    return roots.at(obs)->get_n_evals();
}
