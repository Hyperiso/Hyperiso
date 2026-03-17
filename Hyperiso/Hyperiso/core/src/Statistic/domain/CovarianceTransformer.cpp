#include "CovarianceTransformer.h"

std::vector<std::vector<double>> CovarianceTransformer::transform(const std::vector<ParamId>& ids) {
    std::vector<std::vector<double>> out;

    for (const auto& elem : ids) {
        std::vector<double> temp;
        for (const auto& elem2 : ids) {
            double rho = (*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED);
            temp.push_back(rho);
        }
        out.push_back(std::move(temp));
    }

    return out;
}

std::map<ParamId, std::map<ParamId, double>>
CovarianceTransformer::transform(const std::map<ParamId, double>& ids) {
    std::map<ParamId, std::map<ParamId, double>> out;

    for (const auto& elem : ids) {
        for (const auto& elem2 : ids) {
            out[elem.first][elem2.first] =
                (*corr_proxy)(elem.first, elem2.first,
                              CorrelationProvider::CorrelationType::COMBINED); // TODO Theo : correlation for copula
        }
    }

    return out;
}

std::vector<std::vector<double>>
CovarianceTransformer::transform(const std::vector<ExperimentObs>& ids) {
    std::vector<std::vector<double>> out;

    for (const auto& elem : ids) {
        std::vector<double> temp;
        for (const auto& elem2 : ids) {
            double rho = (*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED);
            temp.push_back(rho);
        }
        out.push_back(std::move(temp));
    }

    return out;
}

std::map<ExperimentObs, std::map<ExperimentObs, double>>
CovarianceTransformer::transform(const std::map<ExperimentObs, double>& ids) {
    std::map<ExperimentObs, std::map<ExperimentObs, double>> out;

    for (const auto& elem : ids) {
        for (const auto& elem2 : ids) {
            out[elem.first][elem2.first] =
                (*corr_proxy)(elem.first, elem2.first,
                              CorrelationProvider::CorrelationType::COMBINED); // TODO Theo : correlation for copula
        }
    }

    return out;
}

std::map<BinnedObservableId, std::map<BinnedObservableId, double>>
CovarianceTransformer::transform(const std::string& experiment,
                                 const std::map<BinnedObservableId, double>& ids) {
    std::map<BinnedObservableId, std::map<BinnedObservableId, double>> out;

    for (const auto& elem : ids) {
        for (const auto& elem2 : ids) {
            out[elem.first][elem2.first] =
                (*corr_proxy)(experiment, elem.first, elem2.first,
                              CorrelationProvider::CorrelationType::COMBINED); // TODO Theo : correlation for copula
        }
    }

    return out;
}

std::vector<std::vector<double>>
CovarianceTransformer::transform(const std::string& experiment,
                                 const std::vector<BinnedObservableId>& ids) {
    std::vector<std::vector<double>> out;

    for (const auto& elem : ids) {
        std::vector<double> temp;
        for (const auto& elem2 : ids) {
            double rho = (*corr_proxy)(experiment, elem, elem2,
                                       CorrelationProvider::CorrelationType::COMBINED);
            temp.push_back(rho);
        }
        out.push_back(std::move(temp));
    }

    return out;
}

std::vector<ParamId> CovarianceTransformer::check_if_corr(const std::vector<ParamId>& ids) {
    std::vector<ParamId> out;

    for (const auto& elem : ids) {
        if ((*par_proxy)(elem, DataType::STD_COMBINED) > 0) { // TODO : eps 1e-50 ? ?
            out.push_back(elem);
        }
    }

    return out;
}