#include "Compound.h"

// TODO : Ensure compatibility with complex compounds (nonlinear operations).

scalar_t Compound::compute_pdv(const ParamId &param_id) const {
    LOG_DEBUG("Computing pdv wrt", param_id);
    ObsParameterProxy opp = ObsParameterProxy();
    scalar_t h = opp(param_id) * 1e-5;
    h = fpeq(std::abs(h), 0.) ? scalar_t(1e-5) : h;

    ObsParameterMutator opm;
    opm.change_mode(param_id, ParameterMode::SHIFTABLE);
    opm.shift(param_id, h);
    scalar_t f_p = eval();
    opm.shift(param_id, -h);
    opm.change_mode(param_id, ParameterMode::FIXED);
    return (f_p - central_value) / h;
}

void Compound::update_gradient() {
    LOG_DEBUG("Updating gradient");
    central_value = eval();
    for (auto &&p : dependences) {
        gradient.insert_or_assign(p, compute_pdv(p));
    }
}

std::unordered_set<ParamId> Compound::get_common_dependences_with(const Compound &other) const {
    std::unordered_set<ParamId> common;
    std::set_intersection(dependences.begin(), dependences.end(),
                          other.get_dependences().begin(), other.get_dependences().end(),
                          std::inserter(common, common.begin()));
    return common;
}

void Compound::add_dependence(const ParamId &param_name) {
    dependences.emplace(param_name);
    if (central_value == NAN) {
        central_value = eval();
    }
    gradient.emplace(param_name, compute_pdv(param_name));
}

void Compound::add_dependences(const std::unordered_set<ParamId> &param_names) {
    LOG_DEBUG("Adding parameter list to compound");
    std::set_union(dependences.begin(), dependences.end(),
                   param_names.begin(), param_names.end(),
                   std::inserter(dependences, dependences.begin()));
    update_gradient();
}

const std::unordered_set<ParamId> &Compound::get_dependences() const {
    return dependences;
}

const std::unordered_map<ParamId, scalar_t> &Compound::get_gradient() const {
    return gradient;
}

scalar_t Compound::variance() {
    scalar_t var = 0;
    ObsParameterProxy opp = ObsParameterProxy();
    CorrelationProxy cp = CorrelationProxy();
    for (const auto &pid_1 : dependences) {
        for (const auto &pid_2 : dependences) {
            if (pid_1 == pid_2) {
                var += pow(opp(pid_1, DataType::STD_COMBINED) * gradient.at(pid_1), 2);
            } else {
                var += cp(pid_1, pid_2, CorrelationProvider::CorrelationType::COMBINED)  // rho_12
                        * opp(pid_1, DataType::STD_COMBINED)                          // sigma_1
                        * opp(pid_2, DataType::STD_COMBINED)                          // sigma_2
                        * gradient.at(pid_1)                                                             // dC/dp_1   
                        * gradient.at(pid_2);                                                            // dC/dp_2
            }
        }
    }
    LOG_DEBUG("Computing compound variance =", var);
    return var;
}

const std::unordered_map<ParamId, scalar_t> Compound::get_leading_uncertainties(size_t n) const {
    std::unordered_map<ParamId, scalar_t> uncertainties = get_uncertainties();
    std::unordered_map<ParamId, scalar_t> max_uncertainties;
    for (size_t i = 0; i < n; i++) {
        std::pair<ParamId, scalar_t> max {{ParameterType::SM, "", 0}, 0};
        for (auto &&p : uncertainties) {
            if (!max_uncertainties.contains(p.first) && p.second > max.second) {
                max = std::make_pair(p.first, p.second);
            }
        }
        max_uncertainties.emplace(max);
    }

    return max_uncertainties;
}

const std::unordered_map<ParamId, scalar_t> Compound::get_uncertainties() const {
    std::unordered_map<ParamId, scalar_t> uncertainties;
    ObsParameterProxy opp = ObsParameterProxy();
    for (auto p : dependences) {
        scalar_t u = opp(p, DataType::STD_COMBINED) * std::abs(gradient.at(p));
        uncertainties.emplace(p, u);
    }

    return uncertainties;
}

Estimate Compound::get_estimate() {
    // TODO : Manage stat and syst separation
    return Estimate {this->eval(), std::sqrt(this->variance()), 0};
}

scalar_t Compound::correlation_with(const Compound &other) const {
    scalar_t corr = 0;
    auto common_dep = get_common_dependences_with(other);
    CorrelationProxy cprox;
    ObsParameterProxy opp = ObsParameterProxy();
    CorrelationProxy cp = CorrelationProxy();
    for (auto &&p_1 : common_dep) {
        for (auto &&p_2 : common_dep) {
            corr += cp(p_1, p_2, CorrelationProvider::CorrelationType::COMBINED)   // rho_12
                    * opp(p_1, DataType::STD_COMBINED)                          // sigma_1
                    * opp(p_2, DataType::STD_COMBINED)                          // sigma_2
                    * this->gradient.at(p_1)                                                       // dC_1/dp_1   
                    * other.get_gradient().at(p_2);                                                // dC_2/dp_2
        }
    }
    return corr;
}

void Compound::print_gradient(std::ostream& os) const {
    for (auto &&p: gradient) {
        os << "d/d(" << p.first.block << "," << p.first.code << ") = " << p.second << std::endl;
    }
}
