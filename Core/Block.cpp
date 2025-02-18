#include "Block.h"

double MapBlock::getValue(int pdgCode) const {
    auto it = values.find(pdgCode);
    if (it != values.end()) {
        return it->second.get_val();
    }
    throw std::invalid_argument("PDG code not found in " + this->blockname);
}

void MapBlock::setValue(int pdgCode, double value, bool force) {
    JSONParser::getInstance(0)->addElement(this->blockname.substr(0, this->blockname.size()-5), pdgCode, value);
    Parameter param (ParamId {ParameterType::CUSTOM, this->blockname.substr(0, this->blockname.size()-5), pdgCode}, value, 0);
    if (force) {
        values[pdgCode] = param;
    } else {
        values[pdgCode] = param;
    }
}

void MapBlock::setMode(int pdgCode, ParameterMode mode) {
    values.at(pdgCode).set_mode(mode);
}

std::map<int, double> MapBlock::getAllValues() {
    std::map<int, double> map_values;
    for (auto& value : values) {
        map_values[value.first] = value.second.get_val();
    }
    return map_values;
}

double MassBlock::getValue(int pdgCode) const {
    if (pdgCode == 5 || pdgCode == 6) {
        LOG_WARN("Accessing heavy quark masses through Parameters is deprecated. Use QCDHelper instead.");
    }
    return MapBlock::getValue(pdgCode);
}

double WilsonBlock::getValue(int pdgCode) const {
    if (pdgCode == -1) {
        return scale;
    } else if (pdgCode == -2) {
        return type;
    }

    int order = pdgCode % 10; 
    WCoef id = static_cast<WCoef>((pdgCode - order) / 10); 

    return values.at(id)[order];
}

void WilsonBlock::setValue(int pdgCode, double value, bool force) {
    if (pdgCode == -1) {
        scale = value;
    } else if (pdgCode == -2) {
        type = (int)value;
    }

    int order = pdgCode % 10; 
    WCoef id = static_cast<WCoef>((pdgCode - order) / 10); 
    if (!values.contains(id)) {
        values.emplace(std::make_pair(id, std::array<double, 3>()));
    }
    values.at(id)[order] = value;
}