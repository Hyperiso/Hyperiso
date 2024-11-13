#include "ModelEvaluator.h"
#include "json_parser.h"
#include "config.hpp"

ModelEvaluator::ModelEvaluator() {}

ModelEvaluator::ModelEvaluator(const std::vector<std::shared_ptr<Observable>>& observables) : observables(observables) {
    LOG_INFO("Instantiating chi squared calculator with", this->observables.size(), "observables");
    update_exp_covariance();
    LOG_INFO("Experimental covariance updated");
    printMatrix(exp_cov_mtx, getDiagonalElements(exp_cov_mtx));
    update_th_covariance();
    LOG_INFO("Theoretical covariance updated");
    printMatrix(th_cov_mtx, getDiagonalElements(th_cov_mtx));
    
}

// void ModelEvaluator::add_observable(Observable obs) {
//     observables.emplace_back(std::move(obs));
//     update_exp_covariance();
//     update_th_covariance();
// }

// void ModelEvaluator::add_observables(std::vector<Observable> obs) {
//     for (auto &&o : obs)
//         observables.emplace_back(o);

//     update_exp_covariance();
//     update_th_covariance();
// }

std::shared_ptr<Observable> ModelEvaluator::find_from_id(Observables id) {
    for (auto o : observables) {
        if (o->getId() == id)
            return o;
    }
    throw std::runtime_error("Observable not found in list.");
}

void ModelEvaluator::update_th_covariance() {
    for (size_t i = 0; i < observables.size(); i++) {
        for (size_t j = i; j < observables.size(); j++) {
            if (i == j) {
                th_cov_mtx.insert(std::make_pair(std::make_pair(observables.at(i)->getId(), observables.at(i)->getId()), 
                                                 observables.at(i)->variance()));
            } else {
                double corr = observables.at(i)->correlation_with(*observables.at(j));
                th_cov_mtx.insert(std::make_pair(std::make_pair(observables.at(i)->getId(), observables.at(j)->getId()), corr));
                th_cov_mtx.insert(std::make_pair(std::make_pair(observables.at(j)->getId(), observables.at(i)->getId()), corr));
            }
        }
    }
}

void ModelEvaluator::update_exp_covariance() {
    // Fill diagonal elements with exp. variance for each obs
    for (size_t i = 0; i < observables.size(); i++) {
        exp_cov_mtx.insert(std::make_pair(std::make_pair(observables.at(i)->getId(), observables.at(i)->getId()), 
                                          observables.at(i)->get_exp_var()));
    }

    // Read any nonzero correlation between obs pairs from data file and add it to the matrix
    std::vector<Correlation> correlations;
    std::vector<Value> _;
    std::string root = project_root.data();
    read_json(root + "/DataBase/data_exp.json", _, correlations);
    auto mapper = ObservableMapper::GetInstance();
    for (const auto& corr : correlations) {
        auto id_1 = mapper->getObservable(corr.name1);
        auto id_2 = mapper->getObservable(corr.name2);
        try {
            auto o_1 = find_from_id(id_1);
            auto o_2 = find_from_id(id_2);
            auto cov = corr.value * std::sqrt(o_1->get_exp_var() * o_2->get_exp_var());
            exp_cov_mtx.insert(std::make_pair(std::make_pair(o_1->getId(), o_2->getId()), cov));
            exp_cov_mtx.insert(std::make_pair(std::make_pair(o_2->getId(), o_1->getId()), cov));
        } catch (std::runtime_error& e) {
            LOG_WARN("Found correlation between observables ", corr.name1, "and", corr.name2, 
                     ", but they do not contribute to the chi squared estimation.");
        }
    }
}

SparseMatrix<Observables> ModelEvaluator::get_covariance() const {
    SparseMatrix<Observables> cov = th_cov_mtx;
    LOG_INFO("Theoretical covariance matrix size", cov.size());
    auto mapper = ObservableMapper::GetInstance();
    for (auto &&p : exp_cov_mtx) {
        LOG_INFO(mapper->getString(p.first.first), mapper->getString(p.first.second));
        if (cov.contains(p.first)) {
            cov.at(p.first) += p.second;
        } else {
            cov.emplace(p);
        }
    }
    return cov;
}

double ModelEvaluator::chi2() const {

    SparseMatrix<Observables> covariance_mtx = get_covariance();
    SparseMatrix<Observables> precision_mtx = invertMatrix(covariance_mtx, getDiagonalElements(covariance_mtx));

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
