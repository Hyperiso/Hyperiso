#include "Block.h"

void Block::addObserver(std::shared_ptr<Block> observer) {
    observers.push_back(observer);
}

void Block::removeObserver(std::shared_ptr<Block> observer) {
    observers.erase(std::find(observers.begin(), observers.end(), observer));
}

void Block::notifyObservers() {
    for (auto& observer : observers) {
        LOG_DEBUG("Notifying observer", observer->blockname, "from source block", blockname);
        if (observer == nullptr) {
            removeObserver(observer);
            continue;
        }
        observer->update();
    }
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
    if (!this->contains(id)) {
        LOG_ERROR("KeyError", "Block", this->blockname, "doesn't contain parameter", id.to_string(), ". Available keys are", std::make_shared<Block>(*this));
    }
    return this->items.at(id);
}

void Block::store(const LhaID& id, std::shared_ptr<Parameter> param) {
    if (this->contains(id)) {
        LOG_WARN("Block", blockname, "already contains a parameter with id", id.to_string());
    } else {
        this->items.emplace(id, param);
    }
}

void Block::assign(const LhaID& key, std::shared_ptr<Parameter> param) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    *(this->items.at(key)) = *param;
    LOG_DEBUG("Call to notifyObservers from Block::assign(const LhaID&, std::shared_ptr<Parameter>) of", blockname);
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

void Block::store_or_assign(const LhaID &key, std::shared_ptr<Parameter> param) {
    if (this->contains(key)) {
        this->assign(key, param);
    } else {
        this->store(key, param);
    }
}

bool Block::contains(const LhaID& id) const {
    return this->items.contains(id);
}

void Block::remove(const LhaID& id) {
    if (!this->items.contains(id)) {
        LOG_ERROR("Cannot remove non-existing parameter", id, "in block", blockname);
    }

    this->items.erase(id);
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
}

void Block::clear_above() {
    for (auto& param: this->items) {
        param.second->clear_above();
    }
}

void Block::clear_below() {
    for (auto& param: this->items) {
        param.second->clear_below();
    }
}

bool DependentBlock::dependsOn(const std::string& blockName) {
    return sourceBlocks.contains(blockName);
}

void DependentBlock::init() {
    self = shared_from_this();
    if (self) {
        LOG_DEBUG("Adding observer to", sourceBlocks.size(), "source blocks");
        for (auto src : sourceBlocks){
            LOG_DEBUG(src.second->blockname);
            src.second->addObserver(self);   
        }
    } else {
        std::cerr << "Error: DependentBlock must be created with std::make_shared!" << std::endl;
    }
}

void DependentBlock::update() {
    if (frozen) {
        LOG_DEBUG("DependentBlock is frozen, skipping update");
        update_at_unfreeze = true;
    } else if (recalculateLambda 
        && std::all_of(sourceBlocks.begin(), sourceBlocks.end(), 
                       [](std::pair<std::string, std::shared_ptr<Block>> block) { return block.second; })) 
    {
        LOG_DEBUG("Updating dependent block", blockname);
        if (auto self = shared_from_this()) { 
            recalculateLambda(sourceBlocks, self);
        } else {
            std::cerr << "Error: shared_from_this() failed in update()" << std::endl;
        }
    } else {
        LOG_ERROR("Error", "DependentBlock::update() called without all source blocks being set");
    }
    LOG_DEBUG("Call to notifyObservers from DependentBlock::update() of", blockname);
    notifyObservers();
}

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
    LOG_DEBUG("Destruct dependentBlock at", self.get());
    if (self) {
        for (auto src : sourceBlocks){
            src.second->removeObserver(self);   
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

void DependentBlock::assign(const LhaID &key, double value) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    this->items.at(key)->set_expected(value);
}

void DependentBlock::clear_above() {
    for (auto &[name, block] : sourceBlocks) {
        block->removeObserver(self);
    }
}

void DependentBlock::clear_below() {
    this->clear_above();
    std::cout << "clear : " << this->get_name() << std::endl;
    //base case
    if (observers.empty()) {
        return; //no need to destroy itself, if no ref, shared_ptr do it
    }

    for (auto &obs : observers) {
        obs->clear_below();
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