#include "CovarianceTransformer.h"

std::vector<std::vector<double>> CovarianceTransformer::transform(std::vector<ParamId> ids) {
    std::vector<std::vector<double>> out;
    for (auto& elem : ids) {
        std::vector<double> temp;
        for (auto& elem2 : ids) {
            temp.push_back(1);
        }
        out.push_back(temp);
    }
}


std::vector<std::vector<double>> CovarianceTransformer::transform(std::vector<ObservableId> ids) {
    std::vector<std::vector<double>> out;
    for (auto& elem : ids) {
        std::vector<double> temp;
        for (auto& elem2 : ids) {
            temp.push_back(1);
        }
        out.push_back(temp);
    }
}