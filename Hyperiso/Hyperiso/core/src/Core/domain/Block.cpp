#include "Block.h"

void Block::addObserver(std::shared_ptr<Block> observer) {
    observers.push_back(observer);
}

void Block::removeObserver(std::shared_ptr<Block> observer) {
    observers.erase(std::find(observers.begin(), observers.end(), observer));
}

void Block::notifyObservers() {
    for (auto& observer : observers) {
        observer->update();
    }
}

Block::Block(std::shared_ptr<Block> other) {
    this->copy(other);
}

Parameter& Block::retrieve(const LhaID& id) {

    if (!this->contains(id)) {
        LOG_ERROR("KeyError", "Block", this->blockname, "doesn't contain parameter", id.to_string());
    }

    return this->items.at(id);
}

void Block::store(const LhaID& id, Parameter&& param) {
    if (this->contains(id)) {
        LOG_WARN("Block", blockname, "already contains a parameter with id", id.to_string());
    } else {
        this->items.emplace(id, param);
    }
}

std::unordered_set<LhaID> Block::getAllIDs() {
    return get_keys(this->items);
}

void Block::set_owner(ParameterType type) {
    for (auto& [_, param] : items) {
        param.set_owner(type);
    }
}

bool Block::contains(const LhaID& id) const {
    return this->items.contains(id);
}

void Block::update(const LhaID &id, Parameter&& param) {
    if (!this->contains(id)) {
        LOG_ERROR("KeyError", "Cannot update non-existing parameter", id.to_string(), "in block", this->blockname);
    }

    this->items.at(id) = param;
}

void Block::copy(std::shared_ptr<Block> other) {
    this->items = other->getItems();
    this->blockname = other->blockname;
}

void Block::remove(const LhaID& id) {
    if (!this->items.contains(id)) {
        LOG_ERROR("Cannot remove non-existing parameter", id, "in block", blockname);
    }

    this->items.erase(id);
}

// TODO : Remove WilsonBlock

double WilsonBlock::getValue(LhaID pdgCode) const {
    if ((long)pdgCode == -1) {
        return scale;
    } else if ((long)pdgCode == -2) {
        return type;
    }

    int order = pdgCode % 10; 
    WCoef id = static_cast<WCoef>((pdgCode - order) / 10); 

    return values.at(id)[order];
}

void WilsonBlock::setValue(LhaID pdgCode, double value, bool force) {
    if ((long)pdgCode == -1) {
        scale = value;
    } else if ((long)pdgCode == -2) {
        type = (int)value;
    }

    int order = pdgCode % 10; 
    WCoef id = static_cast<WCoef>((pdgCode - order) / 10); 
    if (!values.contains(id)) {
        values.emplace(std::make_pair(id, std::array<double, 3>()));
    }
    values.at(id)[order] = value;
    notifyObservers();
}