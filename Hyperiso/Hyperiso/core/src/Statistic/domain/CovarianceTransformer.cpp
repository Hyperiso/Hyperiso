#include "CovarianceTransformer.h"

std::vector<std::vector<double>> CovarianceTransformer::transform(std::vector<ParamId> ids) {
    std::vector<std::vector<double>> out;
    for (auto& elem : ids) {
        std::vector<double> temp;
        for (auto& elem2 : ids) {
            double _ = (*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED);
            temp.push_back(_);
        }
        out.push_back(temp);
    }
    return out;
}

std::map<ParamId, std::map<ParamId, double>> CovarianceTransformer::transform(std::map<ParamId, double> ids) {
    std::map<ParamId, std::map<ParamId, double>> out;
    for (auto& elem : ids) {
        for (auto& elem2 : ids) {
            out[elem.first][elem2.first] = (*corr_proxy)(elem.first, elem2.first, CorrelationProvider::CorrelationType::COMBINED); // TODO Theo : correlation for copula
        }
    }
    return out;
}

std::map<ObservableId, std::map<ObservableId, double>> CovarianceTransformer::transform(std::map<ObservableId, double> ids) {
    std::map<ObservableId, std::map<ObservableId, double>> out;
    for (auto& elem : ids) {
        for (auto& elem2 : ids) {
            out[elem.first][elem2.first] = (*corr_proxy)(elem.first, elem2.first, CorrelationProvider::CorrelationType::COMBINED); // TODO Theo : correlation for copula
        }
    }
    return out;
}

std::vector<std::vector<double>> CovarianceTransformer::transform(std::vector<BinnedObservableId> ids) {
    std::vector<std::vector<double>> out;
    for (auto& elem : ids) {
        std::vector<double> temp;
        for (auto& elem2 : ids) {
            double _ = (*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED);
            temp.push_back(_) ;
        }
        out.push_back(temp);
    }
    return out;
}

std::vector<ParamId> CovarianceTransformer::check_if_corr(std::vector<ParamId> ids) {
    std::vector<ParamId> out;
    for (auto& elem : ids) {
        if ((*par_proxy)(elem, DataType::STD_COMBINED) > 0) { // TODO : eps 1e-50 ? ? 
            out.push_back(elem);
        }
    }
    return out;
}