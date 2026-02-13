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

std::shared_ptr<Parameter> Block::retrieve(const LhaID& id) {
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
        auto w = self_weak();
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

//     auto& dst = this->items.at(key);

//     dst->overwrite_payload_from(*param);

//     dst->notifyObservers();

//     LOG_DEBUG("Call to notifyObservers from Block::assign(const LhaID&, std::shared_ptr<Parameter>) of", blockname);
//     // notifyObservers();
// }

void Block::assign(const LhaID& key, std::shared_ptr<Parameter> param) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    auto& dst = this->items.at(key);
    dst->overwrite_payload_from(*param);

    // notifier param observers si tu veux (pas obligatoire pour DependentBlock)
    dst->notifyObservers();

    // MAIS surtout notifier block observers
    LOG_DEBUG("Call to notifyObservers from Block::assign(shared_ptr) of", blockname);
    notifyObservers();
}

void Block::assign(const LhaID& key, scalar_t value) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    auto& p = this->items.at(key);
    p->set_expected(value);   // ça notifie déjà les observers du param

    // IMPORTANT: notifier les observers du BLOCK (DependentBlock observe le block)
    LOG_DEBUG("Call to notifyObservers from Block::assign(const LhaID&, double) of", blockname);
    notifyObservers();
}

// void Block::assign(const LhaID &key, scalar_t value) {
//     if (!this->contains(key)) {
//         LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
//     }

//     auto& p = this->items.at(key);

//     p->set_expected(value);

//     // p->notifyObservers();

//     LOG_DEBUG("Call to notifyObservers from Block::assign(const LhaID&, double) of", blockname);
//     // notifyObservers();
// }

void Block::store_or_assign(const LhaID& id, std::shared_ptr<Parameter> param) {
    if (!contains(id)) {
        store(id, std::move(param));
    } else {
        assign(id, std::move(param));
    }
}

bool Block::contains(const LhaID& id) const {

    const_cast<Block*>(this)->ensure_up_to_date();

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
    ensure_up_to_date();
    return get_keys(this->items);
}

const std::map<LhaID, std::shared_ptr<Parameter>>& Block::getItems() {
    ensure_up_to_date();
    return this->items;
};

void Block::copy(std::shared_ptr<Block> other) {
    this->items = other->getItems();
    this->blockname = other->blockname;
    if (other->has_scale()) this->set_scale(other->get_scale());

    auto w = self_weak();
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

void DependentBlock::freeze() {
    this->frozen = true;
    Block::freeze();
}


void DependentBlock::unfreeze() {
    this->frozen = false;

    Block::unfreeze();

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


void DependentBlock::assign(const LhaID& key, std::shared_ptr<Parameter> param)
{
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    auto& dst = this->items.at(key);
    dst->overwrite_payload_from(*param);
}




void DependentBlock::assign(const LhaID &key, double value) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }
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

void DependentBlock::mark_dirty() {
    if (dirty) return;
    dirty = true;

    for (auto& obs : observers) {
        if (auto dep = std::dynamic_pointer_cast<DependentBlock>(obs)) {
            dep->mark_dirty();
        } else if (obs) {

        }
    }
}

// void DependentBlock::update() {
//     if (frozen) { update_at_unfreeze = true; return; }
//     mark_dirty();
// }

void DependentBlock::update() {
    if (frozen) { update_at_unfreeze = true; return; }

    mark_dirty();

    for (auto& [_, p] : items) {
        if (!p) continue;
        p->notifyObservers(); 
    }
}

void DependentBlock::ensure_up_to_date_impl() {
    if (frozen) return;
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

    for (auto& [_, p] : items) {
        if (p) p->notifyObservers();
    }
}

std::ostream &operator<<(std::ostream &os, std::shared_ptr<Block> ba) {
    os << "Block " << ba->get_name() << ":\n";
    for (auto &[id, val] : ba->getItems()) {
        os << '\t' << id << ": " << val->get_val() << '\n';
    }
    os << '\n';
    return os;
}