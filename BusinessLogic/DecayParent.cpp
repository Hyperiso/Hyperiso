#include "DecayParent.h"

std::shared_ptr<WilsonInterface> DecayParent::get_wilsons(bool force_update) {
    if (!wilson || force_update) {
        wilson = std::make_shared<WilsonInterface>();
        wilson->build(winfo.wgroups, winfo.matching_scale, winfo.hadronic_scale, winfo.order);
        
        if (winfo.basis == BWilsonBasis::TRADITIONAL 
            && std::find(winfo.wgroups.begin(), winfo.wgroups.end(), WGroup::B) != winfo.wgroups.end()) {
            wilson->switchbasis(WGroup::B);
        }
    } else {
        // for (auto &p : manager->getGroups()) {
            // if (p.second->cacheChanged()) {
                // manager->update(p.first);
            // }
        // }
    }

    return wilson;
}

scalar_t DecayParent::compute_observable(Observables obs) {
    roots.at(obs);
    return roots.at(obs)->calculate();
}
