#include "Observable.h"

Observables Observable::getId() const {
    return id;
}

double Observable::get_exp_val() const {
    return exp_val;
}

double Observable::get_exp_var() const {
    return std::pow(exp_std, 2);
}

double Observable::eval() const {
    return NAN;
}

CoefficientManager *Observable::computeWilsons(bool traditional_basis) const {
    return computeWilsons(model, order, scale, traditional_basis);
}

CoefficientManager *Observable::computeWilsons(int model,
                                               int order,
                                               double scale,
                                               bool traditional_basis) const {
    CoefficientManager* manager;
    double m_W = (*Parameters::GetInstance(0))("MASS", 24);

    switch (model) {
        case 0:
            manager = CoefficientManager::Builder(std::to_string(model), std::map<std::string,std::shared_ptr<CoefficientGroup>>({std::make_pair("BCoefficient", std::make_shared<BCoefficientGroup>(81.)),
            std::make_pair("BScalarCoefficient", std::make_shared<BScalarCoefficientGroup>(81.)),std::make_pair("BPrimeCoefficient", std::make_shared<BPrimeCoefficientGroup>(81.))}), m_W, scale, "LO");
            break;
        case 1:
            manager = CoefficientManager::Builder(std::to_string(model), std::map<std::string,std::shared_ptr<CoefficientGroup>>({std::make_pair("BCoefficient", std::make_shared<BCoefficientGroup_susy>(81.)),
            std::make_pair("BCoefficient", std::make_shared<BScalarCoefficientGroup_susy>(81.)),std::make_pair("BPrimeCoefficient", std::make_shared<BPrimeCoefficientGroup_susy>(81.))}), m_W, scale, "LO");
            break;
        case 2:
            manager = CoefficientManager::Builder(std::to_string(model), std::map<std::string,std::shared_ptr<CoefficientGroup>>({std::make_pair("BCoefficient", std::make_shared<BCoefficientGroup_THDM>(81.)),
            std::make_pair("BScalarCoefficient", std::make_shared<BScalarCoefficientGroup_THDM>(81.)),std::make_pair("BPrimeCoefficient", std::make_shared<BPrimeCoefficientGroup_THDM>(81.))}), m_W, scale, "LO");
            break;
        default:
            LOG_ERROR("ModelError", "Unknown model requested for Wilson coefficient calculation.");
            // return nullptr;
    }
    return manager;
}
