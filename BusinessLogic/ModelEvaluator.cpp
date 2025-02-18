#include "ModelEvaluator.h"
#include "json_parser.h"
#include "config.hpp"

ModelEvaluator::ModelEvaluator() {}

ModelEvaluator::ModelEvaluator(const std::vector<std::shared_ptr<Observable>>& observables) : observables(std::move(observables)) {
    LOG_INFO("Instantiating chi squared calculator with", this->observables.size(), "observables");
    update_exp_data();
    LOG_INFO("Experimental covariance updated");
    printMatrix(exp_cov_mtx, getDiagonalElements(exp_cov_mtx));
    update_th_covariance();
    LOG_INFO("Theoretical covariance updated");
    printMatrix(th_cov_mtx, getDiagonalElements(th_cov_mtx));
}

bool ModelEvaluator::has_observable(Observables id) {
    return static_cast<bool>(find_from_id(id));
}

void ModelEvaluator::add_observable(std::shared_ptr<Observable> obs) {
    if (!has_observable(obs->getId())) {
        observables.emplace_back(obs);
    } else {
        LOG_WARN("ModelEvaluator already takes observable", ObservableMapper::str(obs->getId()), "into account.");
    }
}

void ModelEvaluator::remove_observable(Observables id) {
    if (has_observable(id)) {
        observables.erase(std::find(observables.begin(), observables.end(), find_from_id(id)));
    } else {
        LOG_WARN("ModelEvaluator doesn't take observable", ObservableMapper::str(id), "into account.");
    }
}

std::shared_ptr<Observable> ModelEvaluator::find_from_id(Observables id) {
    for (auto o : observables) {
        if (o->getId() == id)
            return o;
    }
    return nullptr;
}

void ModelEvaluator::update_th_covariance() {
    for (size_t i = 0; i < observables.size(); i++) {
        for (size_t j = i; j < observables.size(); j++) {
            auto id_i = observables.at(i)->getId();
            auto id_j = observables.at(j)->getId();
            if (i == j) {
                th_cov_mtx.insert_or_assign(std::make_pair(id_i, id_j), observables.at(i)->variance());
            } else {
                double corr = observables.at(i)->correlation_with(*observables.at(j));
                th_cov_mtx.insert_or_assign(std::make_pair(id_i, id_j), corr);
                th_cov_mtx.insert_or_assign(std::make_pair(id_j, id_i), corr);
            }
        }
    }
}

void ModelEvaluator::update_exp_data() {
    std::vector<Correlation> correlations;
    std::vector<Value> values;
    std::string exp_path = MemoryManager::GetInstance()->getObservableCovariancePath().string();
    read_json(exp_path, values, correlations);
    LOG_DEBUG("JSON input file read");

    // Read central values for the stored observable and fill diagonal elements of cov matrix with exp. variance 
    std::map<Observables, double> stds; 
    for (const auto& val : values) {
        LOG_DEBUG("Found observable", val.name);

        Observables id = ObservableMapper::enum_elt(val.name);
        auto obs = find_from_id(id);
        if (obs) {
            obs->set_exp_val(val.central_value);
            double std = std::hypot(val.stat_error, val.syst_error);
            stds.emplace(id, std);
            exp_cov_mtx.insert(std::make_pair(std::make_pair(id, id), std * std));
        }
    }

    for (const auto& o : observables) {
        if (fpeq(o->get_exp_val(), 0.)) {
            LOG_ERROR("Experimental data not found for observable ", ObservableMapper::str(o->getId()));
        }
    }

    // Read any nonzero correlation between obs pairs from data file and add it to the matrix
    for (const auto& corr : correlations) {
        auto id_1 = ObservableMapper::enum_elt(corr.name1);
        auto id_2 = ObservableMapper::enum_elt(corr.name2);
        auto o_1 = find_from_id(id_1);
        auto o_2 = find_from_id(id_2);
        if (o_1 && o_2) {
            auto cov = corr.value * stds.at(id_1) * stds.at(id_2);
            exp_cov_mtx.insert(std::make_pair(std::make_pair(id_1, id_2), cov));
            exp_cov_mtx.insert(std::make_pair(std::make_pair(id_2, id_1), cov));
        }
    }
}

SparseMatrix<Observables> ModelEvaluator::get_covariance() {
    update_exp_data();
    LOG_DEBUG("Experimental values retrieved");
    update_th_covariance();
    LOG_DEBUG("Theoretical covariance matrix calculated");

    SparseMatrix<Observables> cov = th_cov_mtx;
    LOG_DEBUG("Theoretical covariance matrix size", cov.size());
    for (auto &&p : exp_cov_mtx) {
        if (cov.contains(p.first)) {
            cov.at(p.first) += p.second;
        } else {
            cov.emplace(p);
        }
    }
    return cov;
}

double ModelEvaluator::chi2() {
    SparseMatrix<Observables> covariance_mtx = get_covariance();
    LOG_DEBUG("Covariance calculated");
    SparseMatrix<Observables> precision_mtx = invertMatrix(covariance_mtx, getDiagonalElements(covariance_mtx));
    LOG_DEBUG("Precision matrix calculated");

    double chi2 {0};
    for (auto &&o_i : observables) {
        for (auto &&o_j : observables) {
            double err_i = o_i->eval() - o_i->get_exp_val();
            double err_j = o_j->eval() - o_j->get_exp_val();
            std::pair<Observables, Observables> obs_pair = std::make_pair(o_i->getId(), o_j->getId());
            if (precision_mtx.contains(obs_pair)) {
                double C_ij_inv = precision_mtx.at(obs_pair);
                chi2 += err_i * C_ij_inv * err_j;
            }
        }
    }

    return chi2;
}
