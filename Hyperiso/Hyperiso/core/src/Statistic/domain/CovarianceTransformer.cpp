#include "CovarianceTransformer.h"

std::vector<std::vector<double>> CovarianceTransformer::transform(std::vector<ParamId> ids) {
    std::vector<std::vector<double>> out;
    for (auto& elem : ids) {
        std::vector<double> temp;
        for (auto& elem2 : ids) {
            double _ = (*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED) * (*par_proxy)(elem, DataType::STD_COMBINED) * (*par_proxy)(elem2, DataType::STD_COMBINED);
            temp.push_back(_);
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
            double _ = (*corr_proxy)(elem, elem2, CorrelationProvider::CorrelationType::COMBINED) * (*par_proxy)(elem, DataType::STD_COMBINED) * (*par_proxy)(elem2, DataType::STD_COMBINED);
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