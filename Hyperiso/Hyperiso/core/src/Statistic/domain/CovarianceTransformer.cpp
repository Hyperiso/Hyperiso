#include "CovarianceTransformer.h"

std::vector<std::vector<double>> CovarianceTransformer::transform(std::vector<ParamId> ids) {
    std::vector<std::vector<double>> out;
    for (auto& elem : ids) {
        std::vector<double> temp;
        for (auto& elem2 : ids) {
            temp.push_back((*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED));
        }
        out.push_back(temp);
    }
    return out;
}


std::vector<std::vector<double>> CovarianceTransformer::transform(std::vector<ObservableId> ids) {
    std::vector<std::vector<double>> out;
    for (auto& elem : ids) {
        std::vector<double> temp;
        for (auto& elem2 : ids) {
            temp.push_back((*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED));
        }
        out.push_back(temp);
    }
    return out;
}