#include "Block.h"

void Block::addObserver(std::shared_ptr<Block> observer) {
    observers.push_back(observer);
}

void Block::notifyObservers() {
    for (auto& observer : observers) {
        observer->update();
    }
}

MapBlock::MapBlock(std::shared_ptr<Block> other) {
    this->copy(other);
}

double MapBlock::getValue(LhaID id) const
{
    auto it = values.find(id);
    if (it != values.end()) {
        return it->second.get_val();
    }
    throw std::invalid_argument("PDG code not found in " + this->blockname);
}

void MapBlock::setValue(LhaID id, double value, bool force) {
    // JSONParser::getInstance(0)->addElement(this->blockname.substr(0, this->blockname.size()-5), id, value);
    Parameter param (ParamId {ParameterType::CUSTOM, this->blockname.substr(0, this->blockname.size()-5), id}, value, 0, 0);
    if (force) {
        values[id] = param;
    } else {
        values[id] = param;
    }
    notifyObservers();
}

void MapBlock::setDeviation(LhaID id, double std_stat, double std_syst, bool force) {
    // JSONParser::getInstance(0)->addElement(this->blockname.substr(0, this->blockname.size()-5), id, std_stat);
    values[id].set_std(std_stat, std_syst);
    notifyObservers();
}

void MapBlock::setMode(LhaID id, ParameterMode mode) {
    values.at(id).set_mode(mode);
    notifyObservers();
}

std::map<LhaID, double> MapBlock::getAllValues() {
    std::map<LhaID, double> map_values;
    for (auto& value : values) {
        map_values[value.first] = value.second.get_val();
    }
    return map_values;
}

std::vector<LhaID> MapBlock::getAllIDs() {
    std::vector<LhaID> ids;
    for (auto& [k, _] : values) 
        ids.emplace_back(k);
    return ids;
}

bool MapBlock::hasID(LhaID id) {
    return values.contains(id);
}

void MapBlock::copy(std::shared_ptr<Block> other) {
    this->values = other->getItems();
}

Parameter MapBlock::getParameter(LhaID id) const {
    if (!values.contains(id)) {
        LOG_ERROR("MapBlock", "Unknown parameter", id, "in block", blockname);
    }
    return this->values.at(id);
}

void MapBlock::setParameter(LhaID id, const Parameter &source) {
    this->values.insert_or_assign(id, source);
    notifyObservers();
}

void MapBlock::remove_parameter(LhaID id) {
    if (!values.contains(id)) {
        LOG_ERROR("Cannot remove non-existing parameter", id, "in block", blockname);
    }
    this->values.erase(id);
}

// void DependentBlock::init() {
//     if (auto self = shared_from_this()) {
//         sourceBlock->addObserver(self);
//     } else {
//         std::cerr << "Error: DependentBlock must be created with std::make_shared!" << std::endl;
//     }
// }

// void DependentBlock::update() {
//     if (sourceBlock) {
//         recalculate();
//     } else {
//         this->values.clear();
//     }
// }

// void DependentBlock::recalculate() {
//     std::cout << "Updating dependent block: " << blockname << " based on " 
//               << sourceBlock->blockname << std::endl;

//     if (sourceBlock) {
//         double newVal = sourceBlock->getValue(1) * 2;
//         setValue(1, newVal, true);
//     }
// }

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

// void GaugeBlock::recalculate() {
//     if (sourceBlock) {
//         double e_em = std::sqrt(4 * M_PI / sourceBlock->getValue(1));
//         double theta_W = 0.5 * std::asin(e_em / (std::sqrt(M_SQRT2 * sourceBlock->getValue(2)) * sourceBlock->getValue(4)));
//         setValue(1, e_em / std::sin(theta_W), true);
//         setValue(2, e_em / std::cos(theta_W), true);
//         setValue(3, std::sqrt(4 * M_PI * sourceBlock->getValue(3)), true);
//         setValue(4, e_em, true);
//     }
// }

// void ReCKMBlock::recalculate() {
//     if (sourceBlock) {
//         double lambda   = sourceBlock->getValue(1);
//         double A        = sourceBlock->getValue(2);
//         double rho      = sourceBlock->getValue(3);
//         double l2 = lambda * lambda;
//         double l3 = lambda * l2;
//         setValue({0, 0}, 1 - l2 / 2, true);
//         setValue({0, 1}, lambda, true);
//         setValue({0, 2}, A * l3 * rho, true);
//         setValue({1, 0}, -lambda, true);
//         setValue({1, 1}, 1 - l2 / 2, true);
//         setValue({1, 2}, A * l2, true);
//         setValue({2, 0}, A * l3 * (1 - rho), true);
//         setValue({2, 1}, -A * l2, true);
//         setValue({2, 2}, 1, true);

//         double dlambda = sourceBlock->getParameter(1).get_std();
//         double dA = sourceBlock->getParameter(2).get_std();
//         double drho = sourceBlock->getParameter(3).get_std();
//         double d1 = lambda * dlambda;
//         double d2 = l2 * dA + 2 * A * lambda * dlambda;
//         double d3 = rho * l3 * dA + 3 * l2 * A * rho * dlambda + A * l3 * drho;
//         setDeviation({0, 0}, d1, 0, true);
//         setDeviation({0, 1}, dlambda, 0, true);
//         setDeviation({0, 2}, d3, 0, true);
//         setDeviation({1, 0}, dlambda, 0, true);
//         setDeviation({1, 1}, d1, 0, true);
//         setDeviation({1, 2}, d2, 0, true);
//         setDeviation({2, 0}, d3, 0, true);
//         setDeviation({2, 1}, d2, 0, true);
//         setDeviation({2, 2}, 0, 0, true);
//     }
// }

// void ImCKMBlock::recalculate() {
//     if (sourceBlock) {
//         double l3      = std::pow(sourceBlock->getValue(1), 3);
//         double A        = sourceBlock->getValue(2);
//         double eta      = sourceBlock->getValue(4);
//         setValue({0, 0}, 0, true);
//         setValue({0, 1}, 0, true);
//         setValue({0, 2}, -A * l3 * eta, true);
//         setValue({1, 0}, 0, true);
//         setValue({1, 1}, 0, true);
//         setValue({1, 2}, 0, true);
//         setValue({2, 0}, -A * l3 * eta, true);
//         setValue({2, 1}, 0, true);
//         setValue({2, 2}, 0, true);

//         double dlambda = sourceBlock->getParameter(1).get_std();
//         double dA = sourceBlock->getParameter(2).get_std();
//         double deta = sourceBlock->getParameter(4).get_std();
//         double d3 = eta * l3 * dA + 3 * l3 / sourceBlock->getValue(1) * A * eta * dlambda + A * l3 * deta;
//         setDeviation({0, 0}, 0, 0, true);
//         setDeviation({0, 1}, 0, 0, true);
//         setDeviation({0, 2}, d3, 0, true);
//         setDeviation({1, 0}, 0, 0, true);
//         setDeviation({1, 1}, 0, 0, true);
//         setDeviation({1, 2}, 0, 0, true);
//         setDeviation({2, 0}, d3, 0, true);
//         setDeviation({2, 1}, 0, 0, true);
//         setDeviation({2, 2}, 0, 0, true);
//     }
// }
