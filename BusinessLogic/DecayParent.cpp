#include "DecayParent.h"

std::shared_ptr<CoefficientManager> DecayParent::compute_wilsons() {
    auto manager = CoefficientManager::Builder(ModelMapper::str(winfo.model), 
                                                winfo.wgroups, 
                                                winfo.matching_scale, 
                                                winfo.hadronic_scale, 
                                                OrderMapper::str(winfo.order));

    if (winfo.model != Model::SM && winfo.model != Model::CUSTOM) {
        auto manager_sm = CoefficientManager::Builder(ModelMapper::str(Model::SM), 
                                                        winfo.wgroups, 
                                                        winfo.matching_scale, 
                                                        winfo.hadronic_scale, 
                                                        OrderMapper::str(winfo.order));
        manager = CoefficientManager::Concat(manager, manager_sm);
    }
    
    if (winfo.basis == BWilsonBasis::TRADITIONAL 
        && std::find(winfo.wgroups.begin(), winfo.wgroups.end(), WilsonGroups::BCoefficients) != winfo.wgroups.end()) {
        manager->switchbasis(GroupMapper::str(WilsonGroups::BCoefficients));
    }

    return manager;
}

scalar_t DecayParent::compute_observable(Observables obs) {
    return roots.at(obs)->calculate();
}
