#include "DependentParameter.h"
#include "SourcesView.h"

DependentParameter::DependentParameter(
    ParamId pid,
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> sources_in,
    DepParamUpdateFunc recalculateFunc)
: Parameter(pid, 0, 0, 0),
  sources_raw(std::move(sources_in)),
  sources(std::make_unique<ParamSrc>(sources_raw, ParameterTypeMapper::str(pid.type.value()) + "::" + pid.block)),
  recalculateLambda(std::move(recalculateFunc)),
  frozen(false) { dirty = true; }

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

void DependentParameter::mark_dirty() {
    dirty = true;
}

void DependentParameter::update() {
    if (frozen) {
        update_at_unfreeze = true;
        this->notifyObservers();
        return;
    }
    dirty = true;
    notifyObservers();
}

scalar_t DependentParameter::get_val() const {
    const_cast<DependentParameter*>(this)->ensure_up_to_date();
    return expected;
}

void DependentParameter::ensure_up_to_date() {
    if (!dirty) return;
    if (frozen) return;

    for (const auto& [pid, src] : sources->raw()) {
        if (!src) LOG_ERROR("Error", "Null source param for", pid.code);
        (void)src->get_val();
    }

    if (!recalculateLambda) {
        LOG_ERROR("Error", "DependentParameter has no recalculateLambda");
    }

    dirty = false;
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

void DependentParameter::rebind(
    std::unordered_map<ParamId, std::shared_ptr<Parameter>> new_sources,
    DepParamUpdateFunc new_lambda)
{
    clear_above();

    sources_raw = std::move(new_sources);
    sources = std::make_unique<ParamSrc>(sources_raw,
        ParameterTypeMapper::str(id.type.value()) + "::" + id.block);
    recalculateLambda = std::move(new_lambda);

    if (auto me = self.lock()) {
        for (auto& [_, src] : sources->raw()) {
            if (src) src->addObserver(me);
        }
    }

    dirty = true;
    notifyObservers();
}

DependentParameter::~DependentParameter() {
    LOG_DEBUG("Destruct DependentParameter at", self.lock().get());
    if (auto me = self.lock()) {
        for (auto& [_, src] : (*sources).raw()) {
            src->removeObserver(me); 
        }
    }
}