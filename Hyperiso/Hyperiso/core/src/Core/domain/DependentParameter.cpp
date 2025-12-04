#include "DependentParameter.h"
#include "SourcesView.hpp" //avoid loops

DependentParameter::DependentParameter(
    ParamId pid,
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources_in,
    DepParamUpdateFunc recalculateFunc)
: Parameter(pid, 0, 0, 0),
  sources_raw(std::move(sources_in)),
  // si ParamSrc a un ctor (map, who):
  sources(std::make_unique<ParamSrc>(sources_raw, ParameterTypeMapper::str(pid.type.value()) + "::" + pid.block)),
  recalculateLambda(std::move(recalculateFunc)),
  frozen(false) {}

bool DependentParameter::dependsOn(const ParamId& pid) {
    return (*sources).raw().contains(pid);
}


void DependentParameter::clear_above() {
    if (auto me = self.lock()) {
        for (auto& [_, param] : (*sources).raw()) {
            param->removeObserver(me);
        }
    }
}

void DependentParameter::clear_below() {
    auto dependents = std::exchange(observers, {});

    this->clear_above();

    if (auto host = owner_block.lock()) {
        host->erase_local(this->id.code);
    }

    for (auto& obs : dependents) {
        if (obs) obs->clear_below();
    }
}

void DependentParameter::init() {
    auto me_base = shared_from_this();
    auto me_dep  = std::static_pointer_cast<DependentParameter>(me_base);
    self = me_dep; 

    for (auto& [_, src] : (*sources).raw()) {
        src->addObserver(me_base);
    }
}


// void DependentParameter::update() {
//     if (frozen) {
//         LOG_DEBUG("DependentParameter is frozen, skipping update");
//         this->update_at_unfreeze = true;
//     } else if (recalculateLambda 
//         && std::all_of((*sources).raw().begin(), (*sources).raw().end(),
//                        [](const std::pair<ParamId, std::shared_ptr<Parameter>>& p){ return p.second != nullptr; })) 
//     {
//         LOG_DEBUG("Updating dependent parameter value");
//         auto me_dep = std::static_pointer_cast<DependentParameter>(shared_from_this());
//         recalculateLambda(*sources, me_dep); 
//     } else {
//         LOG_ERROR("Error", "DependentParameter::update() called without all source parameters being set");
//     }
// }

void DependentParameter::update() {
    if (frozen) {
        LOG_DEBUG("DependentParameter is frozen, skipping update");
        update_at_unfreeze = true;
        return;
    }

    if (!recalculateLambda) {
        LOG_ERROR("Error", "DependentParameter::update() called without recalculateLambda set");
        return;
    }

    for (const auto& [pid, param_ptr] : sources->raw()) {
        if (!param_ptr) {
            LOG_ERROR("Error",
                      "DependentParameter::update() called with a null source parameter for ",
                      pid.code);
            return;
        }
    }

    LOG_DEBUG("Updating dependent parameter value");
    auto me_dep = std::static_pointer_cast<DependentParameter>(shared_from_this());
    recalculateLambda(*sources, me_dep);
}


void DependentParameter::freeze() {
    this->frozen = true;
}

void DependentParameter::unfreeze() {
    this->frozen = false;
    if (update_at_unfreeze) {
        update();
        update_at_unfreeze = false;
    }
}





DependentParameter::~DependentParameter() {
    LOG_DEBUG("Destruct DependentParameter at", self.lock().get());
    if (auto me = self.lock()) {
        for (auto& [_, src] : (*sources).raw()) {
            src->removeObserver(me); 
        }
    }
}