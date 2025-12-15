#include "Block.h"
#include "SourcesView.h"

void Block::addObserver(std::shared_ptr<Block> observer) {
    observers.push_back(observer);
}


void Block::removeObserver(std::shared_ptr<Block> observer) {
    auto it = std::find_if(observers.begin(), observers.end(),
        [&](const std::shared_ptr<Block>& p){ return p.get() == observer.get(); });
    if (it != observers.end()) observers.erase(it);
}

// void Block::notifyObservers() {
//     auto snapshot = observers;
//     for (auto& observer : snapshot) {
//         if (!observer) continue;
//         LOG_DEBUG("Notifying observer", observer->blockname, "from source block", blockname);
//         observer->update();
//     }
//     observers.erase(std::remove(observers.begin(), observers.end(), nullptr), observers.end());
// }

void Block::notifyObservers() {
    for (size_t i = 0; i < observers.size(); ++i) {
        auto& observer = observers[i];
        if (!observer) continue;
        LOG_DEBUG("Notifying observer", observer->blockname, "from source block", blockname);
        observer->update();
    }
    observers.erase(std::remove(observers.begin(), observers.end(), nullptr), observers.end());
}

std::vector<std::shared_ptr<Block>> Block::getObservers() const {
    return observers;
}

void Block::update() {
    LOG_DEBUG("In Block::update");
    for (auto& [_, param] : this->items) {
        param->update();
    }
}

void Block::freeze() {
    for (auto& [_, param] : this->items) {
        param->freeze();
    }
}

void Block::unfreeze() {
    for (auto& [_, param] : this->items) {
        param->unfreeze();
    }
}

Block::Block(std::shared_ptr<Block> other) {
    this->copy(other);
}

// std::shared_ptr<Parameter> Block::retrieve(const LhaID& id) {
//     if (!this->contains(id)) {
//         LOG_ERROR("KeyError", "Block", this->blockname, "doesn't contain parameter", id.to_string(), ". Available keys are", std::make_shared<Block>(*this));
//     }
//     if (auto dep = std::dynamic_pointer_cast<DependentBlock>(shared_from_this())) {
//         dep->ensure_up_to_date();
//     }
//     return this->items.at(id);
// }

// std::shared_ptr<Parameter> Block::retrieve(const LhaID& id) {
//     if (auto dep = std::dynamic_pointer_cast<DependentBlock>(shared_from_this())) {
//         dep->ensure_up_to_date();
//     }

//     if (!this->contains(id)) {
//         LOG_ERROR("KeyError", "Block", this->blockname, "doesn't contain parameter",
//                   id.to_string(), ". Available keys are", std::make_shared<Block>(*this));
//     }

//     return this->items.at(id);
// }

std::shared_ptr<Parameter> Block::retrieve(const LhaID& id) {
    // ✅ plus de shared_from_this()
    this->ensure_up_to_date();

    if (!this->contains(id)) {
        LOG_ERROR("KeyError", "Block", this->blockname, "doesn't contain parameter",
                  id.to_string(), ". Available keys are", std::make_shared<Block>(*this));
    }
    return this->items.at(id);
}

void Block::store(const LhaID& id, std::shared_ptr<Parameter> param) {
    if (this->contains(id)) {
        LOG_DEBUG("Block", blockname, "already contains a parameter with id", id.to_string());
    } else {
        auto w = weak_from_this();
        if (!w.expired()) param->set_owner_block(w);
        this->items.emplace(id, std::move(param));
    }
}

void Block::erase_local(const LhaID& id) {
    this->items.erase(id);
}

// void Block::assign(const LhaID& key, std::shared_ptr<Parameter> param) {
//     if (!this->contains(key)) {
//         LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
//     }

//     *(this->items.at(key)) = *param;
//     LOG_DEBUG("Call to notifyObservers from Block::assign(const LhaID&, std::shared_ptr<Parameter>) of", blockname);
//     notifyObservers();
// }

void Block::assign(const LhaID& key, std::shared_ptr<Parameter> param) {
    auto it = items.find(key);
    if (it == items.end()) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(),
                  "in block", this->blockname);
    }

    // (optionnel mais safe) : si l’ancien param est dépendant, il s’était abonné à des sources
    if (it->second) {
        it->second->clear_above();
    }

    // rattache au block
    auto w = weak_from_this();
    if (!w.expired()) param->set_owner_block(w);

    // le point crucial : on remplace le pointeur
    it->second = std::move(param);

    LOG_DEBUG("Call to notifyObservers from Block::assign(const LhaID&, shared_ptr<Parameter>) of", blockname);
    notifyObservers();
}

void Block::assign(const LhaID &key, scalar_t value) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    this->items.at(key)->set_expected(value);
    LOG_DEBUG("Call to notifyObservers from Block::assign(const LhaID&, double) of", blockname);
    notifyObservers();
}

// void Block::store_or_assign(const LhaID &key, std::shared_ptr<Parameter> param) {
//     if (this->contains(key)) {
//         this->assign(key, param);
//     } else {
//         this->store(key, param);
//     }
// }

// void Block::store_or_assign(const LhaID& id, std::shared_ptr<Parameter> param) {
//     auto it = items.find(id);
//     if (it == items.end()) {
//         // première fois : on stocke le pointeur
//         store(id, std::move(param));
//         return;
//     }

//     this->assign(id, param);
//     // // ✅ IMPORTANT : on garde le pointeur existant, on met juste à jour les champs
//     // auto& existing = it->second;
//     // if (!existing) {
//     //     existing = std::move(param);
//     //     return;
//     // }

//     // // valeur + erreurs (si tu veux)
//     // existing->set_expected(param->get_val());
//     // existing->set_std(param->get_std().first, param->get_std().second);  // si tu as cette API
// }

void Block::store_or_assign(const LhaID& id, std::shared_ptr<Parameter> param) {
    if (!contains(id)) {
        store(id, std::move(param));
    } else {
        assign(id, std::move(param)); // rebind => DependentParameter conservé
    }
}

bool Block::contains(const LhaID& id) const {
    return this->items.contains(id);
}

void Block::remove(const LhaID& id) {
    if (!this->items.contains(id)) {
        LOG_ERROR("Cannot remove non-existing parameter", id, "in block", blockname);
    }
    this->items.at(id)->clear_below();
    this->items.erase(id);

    auto dependents = std::exchange(observers, {});
    for (auto& obs : dependents) {
        if (obs) obs->destroy();
    }
}

void Block::set_owner(ParameterType type) {
    for (auto& [_, param] : items) {
        param->set_owner(type);
    }
}

std::unordered_set<LhaID> Block::getAllIDs() {
    return get_keys(this->items);
}

void Block::copy(std::shared_ptr<Block> other) {
    this->items = other->getItems();
    this->blockname = other->blockname;
    if (other->has_scale()) this->set_scale(other->get_scale());

    auto w = weak_from_this();
    if (!w.expired()) {
        for (auto& [id, p] : items) if (p) p->set_owner_block(w);
    }
}

void Block::clear_above() {
    for (auto& param: this->items) {
        param.second->clear_above();
    }
}

void Block::clear_below() {
    std::vector<std::shared_ptr<Parameter>> snapshot;
    snapshot.reserve(items.size());
    for (auto& kv : items) snapshot.push_back(kv.second);

    for (auto& p : snapshot) {
        if (p) p->clear_below();
    }
}


bool Block::has_scale() {
    return this->scale.has_value();
}

void Block::set_scale(double scale) {
    if (this->has_scale()) 
        LOG_ERROR("LogicError", "Cannot set scale of a block which already has one.");

    this->scale.emplace(scale);
}

double Block::get_scale() {
    if (!this->has_scale()) 
        LOG_ERROR("LogicError", "Scale for the block ",  this->get_name() ," has not been set, cannot retrieve it.");
    
    return this->scale.value_or(0.0);
}

void Block::destroy() {
    std::cout << "destroying (block) :" << this->blockname << std::endl;

    auto dependents = std::exchange(observers, {});
    
    clear_below();
    items.clear();


    for (auto& obs : dependents) {
        if (obs) obs->destroy();
    }
}

std::unordered_map<std::string, std::shared_ptr<Block>> Block::get_source_blocks() const {
    return {};
}

void DependentBlock::init() {
    auto me_base = shared_from_this();     
    auto me_dep  = std::static_pointer_cast<DependentBlock>(me_base);
    self = me_dep;     

    for (auto& [_, src] : sourceBlocks) src->addObserver(me_base);
}

// void DependentBlock::update() {
//     if (frozen) {
//         LOG_DEBUG("DependentBlock is frozen, skipping update");
//         update_at_unfreeze = true;
//         return;
//     } else if (recalculateLambda 
//         && std::all_of(sourceBlocks.begin(), sourceBlocks.end(),
//                        [](const std::pair<std::string, std::shared_ptr<Block>>& b){ return b.second != nullptr; })) 
//     {
//         LOG_DEBUG("Updating dependent block", blockname);
//         auto me_dep = std::static_pointer_cast<DependentBlock>(shared_from_this());
//         recalculateLambda(BlockSrc(this->sourceBlocks, this->blockname), me_dep);
//     } else {
//         LOG_ERROR("Error", "DependentBlock::update() called without all source blocks being set");
//     }
//     LOG_DEBUG("Call to notifyObservers from DependentBlock::update() of", blockname);
//     notifyObservers();
// }

void DependentBlock::freeze() {
    this->frozen = true;
}

void DependentBlock::unfreeze() {
    this->frozen = false;
    if (this->update_at_unfreeze) {
        this->update();
        this->update_at_unfreeze = false;
    }
}

DependentBlock::~DependentBlock() {
    LOG_DEBUG("Destruct dependentBlock at", self.lock().get());
    if (auto me_dep = self.lock()) {
        auto me_base = std::static_pointer_cast<Block>(me_dep);
        for (auto& [_, src] : sourceBlocks) {
            src->removeObserver(me_base);
        }
    }
}

void DependentBlock::assign(const LhaID &key, std::shared_ptr<Parameter> param)
{
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    *(this->items.at(key)) = *param;
}

// void DependentBlock::assign(const LhaID &key, double value) {
//     if (!this->contains(key)) {
//         LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
//     }

//     this->items.at(key)->set_expected(value);
// }

void DependentBlock::assign(const LhaID &key, double value) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }
    // this->items.at(key)->set_expected_silent(value);
    this->items.at(key)->set_expected(value);
}


std::unordered_map<std::string, std::shared_ptr<Block>> DependentBlock::get_source_blocks() const {
    return this->sourceBlocks;
}

void DependentBlock::clear_above() {
    if (auto me = self.lock()) {
        for (auto& [_, block] : sourceBlocks) {
            block->removeObserver(me);
        }
    }
}

void DependentBlock::clear_below() {
    this->clear_above();

    if (observers.empty()) return;

    auto snapshot = observers;
    for (auto& obs : snapshot) {
        if (obs) obs->clear_below();
    }
}

bool DependentBlock::dependsOn(const std::string& blockName) {
    return sourceBlocks.contains(blockName) && sourceBlocks.at(blockName) != nullptr;
}

void DependentBlock::destroy() {
    LOG_DEBUG("destroying (depblock) :", this->blockname);

    auto dependents = std::exchange(observers, {});

    clear_above();

    for (auto& kv : items) {
        kv.second->clear_below();
    }
    items.clear();

    for (auto& obs : dependents) {
        if (obs) obs->destroy();
    }
}

//TODO ::
void DependentBlock::mark_dirty() {
    if (dirty) return;
    dirty = true;
    // propager "dirty" vers les blocks dépendants SANS recalcul
    for (auto& obs : observers) {
        if (auto dep = std::dynamic_pointer_cast<DependentBlock>(obs)) {
            dep->mark_dirty();
        } else if (obs) {
            // fallback: si c’est un Block non-dependent, tu peux soit l’ignorer,
            // soit garder l’ancien comportement.
        }
    }
}

// remplace update() par: si appelé via notify, on marque dirty
void DependentBlock::update() {
    if (frozen) { update_at_unfreeze = true; return; }
    mark_dirty();
}

// void DependentBlock::ensure_up_to_date() {
//     if (!dirty) return;
//     dirty = false;
//     if (!recalculateLambda) LOG_ERROR("Error", "DependentBlock has no recalculateLambda");
//     auto me_dep = std::static_pointer_cast<DependentBlock>(shared_from_this());
//     recalculateLambda(BlockSrc(this->sourceBlocks, this->blockname), me_dep);
// }

// void DependentBlock::ensure_up_to_date() {
//     if (!dirty) return;

//     for (auto& [name, src] : sourceBlocks) {
//         if (!src) continue;
//         if (auto depSrc = std::dynamic_pointer_cast<DependentBlock>(src)) {
//             depSrc->ensure_up_to_date();
//         }
//     }

//     dirty = false;
//     if (!recalculateLambda) LOG_ERROR("Error", "DependentBlock has no recalculateLambda");
//     auto me_dep = std::static_pointer_cast<DependentBlock>(shared_from_this());
//     recalculateLambda(BlockSrc(this->sourceBlocks, this->blockname), me_dep);
// }

void DependentBlock::ensure_up_to_date_impl() {
    if (!dirty) return;

    for (auto& [name, src] : sourceBlocks) {
        if (!src) continue;
        if (auto depSrc = std::dynamic_pointer_cast<DependentBlock>(src)) {
            depSrc->ensure_up_to_date();
        }
    }

    dirty = false;

    if (!recalculateLambda) LOG_ERROR("Error", "DependentBlock has no recalculateLambda");
    auto me_dep = std::static_pointer_cast<DependentBlock>(shared_from_this());
    recalculateLambda(BlockSrc(this->sourceBlocks, this->blockname), me_dep);

    // ✅ IMPORTANT: propager l’invalidation vers les DependentParameter dépendants
    for (auto& [_, p] : items) {
        if (p) p->notifyObservers(); // eux vont juste dirty=true (rapide)
    }
    // std::cout << this << std::endl;
}

std::ostream &operator<<(std::ostream &os, std::shared_ptr<Block> ba) {
    os << "Block " << ba->get_name() << ":\n";
    for (auto &[id, val] : ba->getItems()) {
        os << '\t' << id << ": " << val->get_val() << '\n';
    }
    os << '\n';
    return os;
}