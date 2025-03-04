#include "Compound.h"

double Compound::compute_pdv(const ParamId &param_id) const {
    LOG_DEBUG("Computing pdv wrt", param_id);
    auto p = Parameters::GetInstance(param_id.type);
    double h = Parameters::Get(param_id) * 1e-5;
    h = fpeq(h, 0.) ? 1e-8 : h;
    p->changeParameterMode(param_id, ParameterMode::SHIFTABLE);
    p->shiftParameter(param_id, h);
    double f_p = eval();
    p->shiftParameter(param_id, -h);
    p->changeParameterMode(param_id, ParameterMode::FIXED);
    return (f_p - central_value) / h;
}

void Compound::update_gradient() {
    LOG_DEBUG("Updating gradient");
    central_value = eval();
    for (auto &&p : dependences) {
        gradient.insert_or_assign(p, compute_pdv(p));
    }
}

std::vector<ParamId> Compound::get_common_dependences_with(const Compound &other) const {
    std::vector<ParamId> common;
    for (auto &&p: dependences) {
        if (std::find(other.get_dependences().begin(), other.get_dependences().end(), p) != other.get_dependences().end())
            common.emplace_back(p);
    }
    return common;
}

void Compound::add_dependence(const ParamId &param_name) {
    dependences.emplace_back(param_name);
    if (central_value == NAN) {
        central_value = eval();
    }
    gradient.emplace(param_name, compute_pdv(param_name));
}

void Compound::add_dependences(const std::vector<ParamId> &param_names) {
    LOG_DEBUG("Adding parameter list to compound");
    for (auto &&p : param_names) {
        dependences.emplace_back(p);
    }
    update_gradient();
}

const std::vector<ParamId> &Compound::get_dependences() const {
    return dependences;
}

const std::map<ParamId, double> &Compound::get_gradient() const {
    return gradient;
}

double Compound::variance() {
    double var = 0;
    CorrelationRepository cr;
    for (const auto &pid_1 : dependences) {
        for (const auto &pid_2 : dependences) {
            if (pid_1 == pid_2) {
                var += std::pow(pid_1.std * gradient.at(pid_1), 2);
            } else {
                var += cr.get_combined_correlation(pid_1, pid_2) * pid_1.std * pid_2.std * gradient.at(pid_1) * gradient.at(pid_2);
            }
        }
    }
    LOG_DEBUG("Computing compound variance =", var);
    return var;
}

const std::map<ParamId, double> Compound::get_leading_uncertainties(size_t n) const {
    std::map<ParamId, double> uncertainties = get_uncertainties();
    std::map<ParamId, double> max_uncertainties;
    for (size_t i = 0; i < n; i++) {
        std::pair<ParamId, double> max {{ParameterType::SM, "", 0}, 0};
        for (auto &&p : uncertainties) {
            if (!max_uncertainties.contains(p.first) && p.second > max.second) {
                max = std::make_pair(p.first, p.second);
            }
        }
        max_uncertainties.emplace(max);
    }

    return max_uncertainties;
}

const std::map<ParamId, double> Compound::get_uncertainties() const {
    std::map<ParamId, double> uncertainties;
    CorrelationRepository cr;
    for (auto p : dependences) {
        double u = std::sqrt(cr.get_combined_correlation(p, p)) * p.std * std::abs(gradient.at(p));
        uncertainties.emplace(p, u);
    }

    return uncertainties;
}

double Compound::correlation_with(const Compound &other) const {
    double corr = 0;
    auto common_dep = get_common_dependences_with(other);
    CorrelationRepository cr;
    for (auto &&p_1 : common_dep) {
        for (auto &&p_2 : common_dep) {
                corr += cr.get_combined_correlation(p_1, p_2) * p_1.std * p_2.std * gradient.at(p_1) * other.get_gradient().at(p_2);
        }
    }
    return corr;
}

void Compound::print_gradient(std::ostream& os) const {
    for (auto &&p: gradient) {
        os << "d/d(" << p.first.block << "," << p.first.code << ") = " << p.second << std::endl;
    }
}
