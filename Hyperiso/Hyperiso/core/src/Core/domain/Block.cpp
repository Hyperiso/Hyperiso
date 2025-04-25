#include "Block.h"

void Block::addObserver(std::shared_ptr<Block> observer) {
    observers.push_back(observer);
}

void Block::removeObserver(std::shared_ptr<Block> observer) {
    observers.erase(std::find(observers.begin(), observers.end(), observer));
}

void Block::notifyObservers() {
    for (auto& observer : observers) {
        LOG_INFO("Notifying observer", observer->blockname, "from source block", blockname);
        observer->update();
    }
}

void Block::update() {
    LOG_INFO("In Block::update");
    for (auto& [_, param] : this->items) {
        param->update();
    }
}

void Block::freeze() {
    LOG_INFO("In Block::freeze");
    for (auto& [_, param] : this->items) {
        param->freeze();
    }
}

void Block::unfreeze() {
    LOG_INFO("In Block::unfreeze");
    for (auto& [_, param] : this->items) {
        param->unfreeze();
    }
}

Block::Block(std::shared_ptr<Block> other) {
    this->copy(other);
}

std::shared_ptr<Parameter> Block::retrieve(const LhaID& id) {
    if (!this->contains(id)) {
        LOG_ERROR("KeyError", "Block", this->blockname, "doesn't contain parameter", id.to_string());
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
    LOG_INFO("Call to notifyObservers from Block::assign(const LhaID&, std::shared_ptr<Parameter>) of", blockname);
    notifyObservers();
}

void Block::assign(const LhaID &key, double value) {
    if (!this->contains(key)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", key.to_string(), "in block", this->blockname);
    }

    this->items.at(key)->set_expected(value);
    LOG_INFO("Call to notifyObservers from Block::assign(const LhaID&, double) of", blockname);
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
