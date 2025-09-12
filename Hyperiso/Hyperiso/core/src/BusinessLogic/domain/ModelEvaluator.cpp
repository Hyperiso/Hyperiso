#include "ModelEvaluator.h"
#include "config.hpp"

ModelEvaluator::ModelEvaluator() {}

ModelEvaluator::ModelEvaluator(const std::unordered_set<std::shared_ptr<Observable>>& observables) {
    for (auto &&obs : observables) {
        this->observables.emplace(obs->getId(), obs);
    }
    LOG_INFO("Instantiating chi squared calculator with", this->observables.size(), "observables");
    printMatrix(exp_cov_mtx, getDiagonalElements(exp_cov_mtx));
    update_th_covariance();
    update_exp_covariance();
    LOG_INFO("Covariance matrices updated");
    printMatrix(th_cov_mtx, getDiagonalElements(th_cov_mtx));
}

bool ModelEvaluator::has_observable(Observables id) {
    return this->observables.contains(id);
}

void ModelEvaluator::add_observable(std::shared_ptr<Observable> obs) {
    if (!has_observable(obs->getId())) {
        observables.emplace(obs->getId(), obs);
    } else {
        LOG_WARN("ModelEvaluator already takes observable", ObservableMapper::str(obs->getId()), "into account.");
    }
}

void ModelEvaluator::remove_observable(Observables id) {
    if (has_observable(id)) {
        observables.erase(id);
    } else {
        LOG_WARN("ModelEvaluator doesn't take observable", ObservableMapper::str(id), "into account.");
    }
}

void ModelEvaluator::update_th_covariance() {
    // TODO : optimisation possible via map et itérateur pour exploiter la symétrie de la matrice de covariance
    for (auto &[id_1, obs_1] : this->observables) {
        for (auto &[id_2, obs_2] : this->observables) {
            auto pair_id = std::make_pair(id_1, id_2);
            if (id_1 == id_2) {
                th_cov_mtx.insert_or_assign(pair_id, obs_1->variance());
            } else {
                double corr = obs_1->correlation_with(*obs_2);
                th_cov_mtx.insert_or_assign(pair_id, corr);
            }
        }
    }
}

void ModelEvaluator::update_exp_covariance() {
    CorrelationProxy cp;
    ObsParameterProxy opp {ParameterType::OBSERVABLE};
    for (auto &[id_1, obs_1] : this->observables) {
        for (auto &[id_2, obs_2] : this->observables) {
            auto pair_id = std::make_pair(id_1, id_2);
            if (id_1 == id_2) {
                scalar_t var = std::pow(opp("FOBS", LhaID(ObservableMapper::str(id_1)), DataType::STD_COMBINED), 2);
                exp_cov_mtx.insert_or_assign(pair_id, var);
            } else {
                double corr = cp(id_1, id_2, CorrelationProvider::CorrelationType::COMBINED);
                scalar_t sigma_1 = opp("FOBS", ObservableMapper::flha_of(ObservableMapper::to_id(id_1)).value(), DataType::STD_COMBINED);
                scalar_t sigma_2 = opp("FOBS", ObservableMapper::flha_of(ObservableMapper::to_id(id_2)).value(), DataType::STD_COMBINED);
                exp_cov_mtx.insert_or_assign(pair_id, corr * sigma_1 * sigma_2);
            }
        }
    }
}

SparseMatrix<Observables> ModelEvaluator::get_covariance() {
    update_th_covariance();
    LOG_DEBUG("Theoretical covariance matrix calculated");

    SparseMatrix<Observables> cov = th_cov_mtx;
    // customPrintMatrix(cov, getDiagonalElements(cov));
    CorrelationProxy corr_prox;
    LOG_DEBUG("Theoretical covariance matrix size", cov.size());

    for (auto &&p : exp_cov_mtx) {
        // std::cout << ObservableMapper::str(p.first.first) << " " << ObservableMapper::str(p.first.second) << ": "<< p.second << std::endl;
        if (cov.contains(p.first)) {
            cov.at(p.first) += p.second;
        } else {
            cov.emplace(p);
        }
    }
    auto final_cov = removeEmptyRowsAndCols(cov, getDiagonalElements(cov));

    // for (auto& elem : getDiagonalElements(cov)) {
    //     std::cout << "before : " << ObservableMapper::str(elem) << std::endl;
    // }
    
    // for (auto& elem : final_cov.second) {
    //     std::cout << "after : " << ObservableMapper::str(elem) << std::endl;
    // }
    // customPrintMatrix(final_cov.first, final_cov.second);
    return final_cov.first;
}

double ModelEvaluator::chi2() {
    SparseMatrix<Observables> covariance_mtx = get_covariance();
    LOG_DEBUG("Covariance calculated");
    SparseMatrix<Observables> precision_mtx = invertMatrix(covariance_mtx, getDiagonalElements(covariance_mtx));
    LOG_DEBUG("Precision matrix calculated");
    // customPrintMatrix(precision_mtx, getDiagonalElements(precision_mtx));

    double chi2 {0};
    for (auto &&[id_i, obs_i] : this->observables) {
        for (auto &&[id_j, obs_j] : this->observables) {
            double err_i = obs_i->eval() - obs_i->get_exp_val();
            double err_j = obs_j->eval() - obs_j->get_exp_val();
            std::pair<Observables, Observables> obs_pair = std::make_pair(id_i, id_j);
            if (precision_mtx.contains(obs_pair)) {
                double C_ij_inv = precision_mtx.at(obs_pair);
                chi2 += err_i * C_ij_inv * err_j;
            }
        }
    }

    return chi2;
}
