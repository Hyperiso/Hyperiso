#include "Parameter.h"
#include <stdexcept>

Parameter::Parameter(ParamId id, scalar_t mean, scalar_t std_stat, scalar_t std_syst) 
    : id(id), expected(mean), deviation_stat(std_stat), deviation_syst(std_syst), shift(0), mode(ParameterMode::FIXED) {}


void Parameter::set_mode(ParameterMode new_mode) {
    if (new_mode != mode) {
        mode = new_mode;
    }
}

void Parameter::set_std(scalar_t stat, scalar_t syst) {
    this->deviation_stat = stat;
    this->deviation_syst = syst;
}

scalar_t Parameter::get_val() const {
    return this->mode == ParameterMode::FIXED ? expected : expected + shift;
}

void Parameter::set_expected(scalar_t val) {
    this->expected = val;
    notifyObservers();
}

void Parameter::set_id(ParamId id) {
    this-> id = id;
}

scalar_t Parameter::get_combined_std() const {
    return std::hypot(static_cast<double>(deviation_stat),
                  static_cast<double>(deviation_syst));
}

std::pair<scalar_t, scalar_t> Parameter::get_std() const {
    return std::make_pair(deviation_stat, deviation_syst);
}

ParamId Parameter::get_id() const {
    return id;
}

void Parameter::set_owner(ParameterType type) {
    id.set_parameter_type(type);
}

void Parameter::set_shift(scalar_t shift) {
    if (mode == ParameterMode::SHIFTABLE) {
        this->shift = shift;
    } else {
        throw std::runtime_error("Cannot change the value of immutable parameter " + id.block + " " + std::to_string(id.code));
    }
}

Parameter& Parameter::operator=(const Parameter& other) {
    this->id = other.id;
    this->set_expected(other.expected);
    this->deviation_stat = other.deviation_stat;
    this->deviation_syst = other.deviation_syst;
    this->mode = other.mode;
    this->shift = other.shift;
    return *this;
}

bool DependentParameter::dependsOn(const ParamId& pid) {
    return sources.contains(pid);
}

void DependentParameter::init() {
    self = shared_from_this();
    if (self) {
        for (auto src : sources){
            src.second->addObserver(self);   
        }
    } else {
        std::cerr << "Error: DependentBlock must be created with std::make_shared!" << std::endl;
    }
}

void DependentParameter::update() {
    if (frozen) {
        LOG_DEBUG("DependentParameter is frozen, skipping update");
        this->update_at_unfreeze = true;
    } else if (recalculateLambda 
        && std::all_of(sources.begin(), sources.end(), 
                    [](std::pair<ParamId, std::shared_ptr<Parameter>> block) { return block.second; })) 
    {
        LOG_DEBUG("Updating dependent parameter value");
        if (auto self = shared_from_this()) { 
            recalculateLambda(sources, self);
        } else {
            std::cerr << "Error: shared_from_this() failed in update()" << std::endl;
        }
    } else {
        LOG_ERROR("Error", "DependentParameter::update() called without all source parameters being set");
    }
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
    LOG_INFO("Destruct DependentParameter at", self.get());
    if (self) {
        for (auto src : sources){
            src.second->removeObserver(self);   
        }
    }
}

Parameter& Parameter::operator+=(const Parameter& other) {
    this->expected       += other.expected;
    this->deviation_stat = std::hypot(other.deviation_stat, this->deviation_stat);
    this->deviation_syst = std::hypot(other.deviation_syst, this->deviation_syst);
    this->shift          += other.shift;
    return *this;
}

Parameter &Parameter::operator*=(const scalar_t &scale) {
    this->expected       *= scale;
    this->deviation_stat *= std::abs(scale);
    this->deviation_syst *= std::abs(scale);
    this->shift          *= scale;
    return *this;
}
