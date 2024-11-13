#include "Compound.h"
#include "Parameter.h"
#include "Parameters.h"
#include "json_parser.h"
#include <algorithm>

void Compound::read_param_covariance() {
    std::vector<Correlation> correlations;
    std::vector<Value> values;
    std::string root = project_root.data();
    read_json(root + "/DataBase/param_cov.json", values, correlations);

    std::map<ParamId, double> stds;

    // Fill diagonal elements with exp. variance for each parameter
    for (const auto& val : values) {
        std::string del = "|";
        auto split = val.name.find(del);
        ParamId id = std::make_pair(val.name.substr(0, split), std::stoi(val.name.substr(split + 1, val.name.size() - split)));
        if (std::find(dependences.begin(), dependences.end(), id) != dependences.end()) {
            param_corr.insert(std::make_pair(std::make_pair(id, id), std::pow(val.stat_error, 2)));
            stds.emplace(std::make_pair(id, val.stat_error));
        }
    }

    // Read any nonzero correlation between obs pairs from data file and add it to the matrix
    for (const auto& corr : correlations) {
        std::string del = "|";
        auto split = corr.name1.find(del);
        ParamId id_1 = std::make_pair(corr.name1.substr(0, split), std::stoi(corr.name1.substr(split + 1, corr.name1.size() - split)));
        split = corr.name2.find(del);
        ParamId id_2 = std::make_pair(corr.name2.substr(0, split), std::stoi(corr.name2.substr(split + 1, corr.name2.size() - split)));
        auto cov = corr.value * stds.at(id_1) * stds.at(id_2);
        param_corr.insert(std::make_pair(std::make_pair(id_1, id_2), cov));
        param_corr.insert(std::make_pair(std::make_pair(id_2, id_1), cov));
    }
}

double Compound::compute_pdv(const ParamId &param_id) const
{
    for (int i = 0; i < N_PARAM_INSTANCES; i++) {
        auto p = Parameters::GetInstance(i);
        try {
            double h = (*p)(param_id.first, param_id.second) * 1e-3;
            p->changeParameterMode(param_id, ParameterMode::SHIFTABLE);
            p->shiftParameter(param_id, h);
            double f_p = eval();
            p->shiftParameter(param_id, -h);
            double f_m = eval();
            p->changeParameterMode(param_id, ParameterMode::SHIFTABLE);
            return (f_p - f_m) / (2 * h);
        } catch(const std::invalid_argument& e) {}
    }
}

void Compound::update_gradient() {
    for (auto &&p : dependences) {
        gradient.emplace(p, compute_pdv(p));
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

void Compound::add_dependence(const ParamId &param_name)
{
    dependences.emplace_back(param_name);
    gradient.emplace(param_name, compute_pdv(param_name));
}

void Compound::add_dependences(const std::vector<ParamId> &param_names) {
    for (auto &&p : param_names) 
        add_dependence(p);
}

const std::vector<ParamId> &Compound::get_dependences() const {
    return dependences;
}

const std::map<ParamId, double> &Compound::get_gradient() const {
    return gradient;
}

double Compound::variance() const {
    double var = 0;
    for (auto &&pp : param_corr) {
        var += pp.second * gradient.at(pp.first.first) * gradient.at(pp.first.second);
    }
    return var;
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
        os << "d/d(" << p.first.first << "," << p.first.second << ") = " << p.second << std::endl;
    }
}
