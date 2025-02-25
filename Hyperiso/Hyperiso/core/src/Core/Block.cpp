#include "Block.h"

double MapBlock::getValue(LhaID id) const {
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
}

void MapBlock::setDeviation(LhaID id, double std_stat, double std_syst, bool force) {
    // JSONParser::getInstance(0)->addElement(this->blockname.substr(0, this->blockname.size()-5), id, std_stat);
    values[id].set_std(std_stat, std_syst);
}

void MapBlock::setMode(LhaID id, ParameterMode mode) {
    values.at(id).set_mode(mode);
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
}