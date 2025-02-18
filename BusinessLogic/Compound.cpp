#include "Compound.h"
#include "Parameter.h"
#include "Parameters.h"
#include "json_parser.h"
#include <algorithm>

void Compound::read_param_covariance() {
    LOG_DEBUG("Extracting parameter covariance from input files.");
    std::vector<Correlation> correlations;
    std::vector<Value> values;
    std::string exp_path = MemoryManager::GetInstance()->getParameterCovariancePath().string();
    read_json(exp_path, values, correlations);

    std::map<ParamId, double> stds;

    // Fill diagonal elements with exp. variance for each parameter
    for (const auto& val : values) {
        std::string del = "|";
        auto split = val.name.find(del);
        std::string block = val.name.substr(0, split);
        int pdg = std::stoi(val.name.substr(split + 1, val.name.size() - split));
        ParameterType type = Parameters::GetType(block, pdg);
        ParamId id {type, block, pdg};
        if (std::find(dependences.begin(), dependences.end(), id) != dependences.end()) {
            param_corr.insert(std::make_pair(std::make_pair(id, id), std::pow(val.stat_error, 2)));
            stds.emplace(std::make_pair(id, val.stat_error));
        }
    }

    // Read any nonzero correlation between parameter pairs from data file and add it to the matrix
    for (const auto& corr : correlations) {
        std::string del = "|";
        auto split = corr.name1.find(del);
        std::string block = corr.name1.substr(0, split);
        int pdg = std::stoi(corr.name1.substr(split + 1, corr.name1.size() - split));
        ParameterType type = Parameters::GetType(block, pdg);
        ParamId id_1 = {type, block, pdg};
        split = corr.name2.find(del);
        block = corr.name2.substr(0, split);
        pdg = std::stoi(corr.name2.substr(split + 1, corr.name2.size() - split));
        type = Parameters::GetType(block, pdg);
        ParamId id_2 = {type, block, pdg};
        auto cov = corr.value * stds.at(id_1) * stds.at(id_2);
        param_corr.insert(std::make_pair(std::make_pair(id_1, id_2), cov));
        param_corr.insert(std::make_pair(std::make_pair(id_2, id_1), cov));
    }
}

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
    read_param_covariance();
}

void Compound::add_dependences(const std::vector<ParamId> &param_names) {
    LOG_DEBUG("Adding parameter list to compound");
    for (auto &&p : param_names) {
        dependences.emplace_back(p);
    }
    update_gradient();
    read_param_covariance();
}

const std::vector<ParamId> &Compound::get_dependences() const {
    return dependences;
}

const std::map<ParamId, double> &Compound::get_gradient() const {
    return gradient;
}

double Compound::variance() {
    double var = 0;
    for (auto &&pp : param_corr) {
        var += pp.second * gradient.at(pp.first.first) * gradient.at(pp.first.second);
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
    for (auto p : dependences) {
        if (param_corr.contains({p, p})) {
            std::pair<ParamId, double> u = {p, std::sqrt(param_corr.at({p, p})) * std::abs(gradient.at(p))};
            uncertainties.emplace(u);
        }
    }

    return uncertainties;
}

double Compound::correlation_with(const Compound &other) const {
    double corr = 0;
    auto common_dep = get_common_dependences_with(other);
    for (auto &&p_1 : common_dep) {
        for (auto &&p_2 : common_dep) {
            if (param_corr.contains({p_1, p_2}))
                corr += param_corr.at({p_1, p_2}) * gradient.at(p_1) * other.get_gradient().at(p_2);
        }
    }
    return corr;
}

void Compound::print_gradient(std::ostream& os) const {
    for (auto &&p: gradient) {
        os << "d/d(" << p.first.block << "," << p.first.code << ") = " << p.second << std::endl;
    }
}
