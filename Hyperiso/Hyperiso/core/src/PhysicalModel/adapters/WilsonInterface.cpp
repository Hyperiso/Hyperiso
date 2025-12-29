#include "WilsonInterface.h"

QCDOrder WilsonInterface::ensure_mty_compat(QCDOrder order) {
    if (UseMarty().get() && !(order == QCDOrder::LO)) {
        LOG_WARN("Using MARTY defaults all calculations to LO in QCD.");
        return QCDOrder::LO;
    }
    return order;
}


void WilsonInterface::build(WilsonBuildConfig config) {
    this->builder = std::make_shared<WilsonBuilder>(config);
    this->provider = this->builder->get_wilson_provider();
    built = true;
}

void WilsonInterface::addWilsonGroup(WilsonBuildConfig config) {
    if (!this->builder) {
        LOG_ERROR("AccessError", "Please build the interface first before adding groups");
    }
    this->builder->add(config);
}

void WilsonInterface::set_matching_scale(double mu_W) {
    ScaleSetter(ScaleType::MATCHING).set(mu_W);
}

void WilsonInterface::set_hadronic_scale(double mu_h) {
    ScaleSetter(ScaleType::HADRONIC).set(mu_h);
}

scalar_t WilsonInterface::getMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
    if (!built) {
        LOG_ERROR("LogicError", "Interface has not been built");
    }

    WilsonRequest request {
        group,
        coeff,
        order,
        cont_type,
        ScaleType::MATCHING,
        false // sum_qcd_orders
    };
    return this->provider->get(std::make_shared<WilsonRequest>(request));
}

scalar_t WilsonInterface::getM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
    return getMatchingCoefficient(group, coeff, order, cont_type);
}

scalar_t WilsonInterface::getFullMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
    if (!built) {
        LOG_ERROR("LogicError", "Interface has not been built");
    }

    WilsonRequest request {
        group,
        coeff,
        order,
        cont_type,
        ScaleType::MATCHING,
        true // sum_qcd_orders
    };
    return this->provider->get(std::make_shared<WilsonRequest>(request));
}

scalar_t WilsonInterface::getFM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type) {
    return getFullMatchingCoefficient(group, coeff, order, cont_type);
}

scalar_t WilsonInterface::getRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    if (!built) {
        LOG_ERROR("LogicError", "Interface has not been built");
    }

    WilsonRequest request {
        group,
        coeff,
        order,
        cont_type,
        ScaleType::HADRONIC,
        false // sum_qcd_orders
    };
    request.basis = basis;
    return this->provider->get(std::make_shared<WilsonRequest>(request));
}

scalar_t WilsonInterface::getR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    return getRunCoefficient(group, coeff, order, cont_type);
}

scalar_t WilsonInterface::getFullRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    if (!built) {
        LOG_ERROR("LogicError", "Interface has not been built");
    }

    WilsonRequest request {
        group,
        coeff,
        order,
        cont_type,
        ScaleType::HADRONIC,
        true // sum_qcd_orders
    };
    request.basis = basis;
    return this->provider->get(std::make_shared<WilsonRequest>(request));
}

scalar_t WilsonInterface::getFR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    return getFullRunCoefficient(group, coeff, order, cont_type);
}

std::map<QCDOrder, scalar_t> WilsonInterface::getSepOrderMatchingCoefficient(WGroup group, WCoef coeff, ContributionType cont_type) {
    std::map<QCDOrder, scalar_t> C {{
        {QCDOrder::LO, getMatchingCoefficient(group, coeff, QCDOrder::LO, cont_type)},
        {QCDOrder::NLO, getMatchingCoefficient(group, coeff, QCDOrder::NLO, cont_type)},
        {QCDOrder::NNLO, getMatchingCoefficient(group, coeff, QCDOrder::NNLO, cont_type)}
    }};
    return C;
}

std::map<QCDOrder, scalar_t> WilsonInterface::getSM(WGroup group, WCoef coeff, ContributionType cont_type) {
    return getSepOrderMatchingCoefficient(group, coeff, cont_type);
}

std::map<QCDOrder, scalar_t> WilsonInterface::getSepOrderRunCoefficient(WGroup group, WCoef coeff, ContributionType cont_type, WilsonBasis basis) {
    std::map<QCDOrder, scalar_t> C {{
        {QCDOrder::LO, getRunCoefficient(group, coeff, QCDOrder::LO, cont_type)},
        {QCDOrder::NLO, getRunCoefficient(group, coeff, QCDOrder::NLO, cont_type)},
        {QCDOrder::NNLO, getRunCoefficient(group, coeff, QCDOrder::NNLO, cont_type)}
    }};
    return C;
}

std::map<QCDOrder, scalar_t> WilsonInterface::getSR(WGroup group, WCoef coeff, ContributionType cont_type, WilsonBasis basis) {
    return getSepOrderRunCoefficient(group, coeff, cont_type);
}

std::map<WCoef, scalar_t> WilsonInterface::getAllMatchingCoefficients(WGroup group, QCDOrder order, ContributionType cont_type) {
    std::vector<WCoef> ids = WCoefMapper::get_group(group);
    std::map<WCoef, scalar_t> Cs;
    for (auto c : ids) {
        Cs.emplace(c, getMatchingCoefficient(group, c, order, cont_type));
    }

    return Cs;
}

std::map<WCoef, scalar_t> WilsonInterface::getAM(WGroup group, QCDOrder order, ContributionType cont_type) {
    return getAllMatchingCoefficients(group, order, cont_type);
}

std::map<WCoef, scalar_t> WilsonInterface::getAllRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    std::vector<WCoef> ids = WCoefMapper::get_group(group);
    std::map<WCoef, scalar_t> Cs;
    for (auto c : ids) {
        Cs.emplace(c, getRunCoefficient(group, c, order, cont_type));
    }

    return Cs;
}

std::map<WCoef, scalar_t> WilsonInterface::getAR(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    return getAllRunCoefficients(group, order, cont_type);
}

std::map<WCoef, scalar_t> WilsonInterface::getAllFullMatchingCoefficients(WGroup group, QCDOrder order, ContributionType cont_type) {
    std::vector<WCoef> ids = WCoefMapper::get_group(group);
    std::map<WCoef, scalar_t> Cs;
    for (auto c : ids) {
        Cs.emplace(c, getFullMatchingCoefficient(group, c, order, cont_type));
    }

    return Cs;
}

std::map<WCoef, scalar_t> WilsonInterface::getAFM(WGroup group, QCDOrder order, ContributionType cont_type) {
    return getAllFullMatchingCoefficients(group, order, cont_type);
}

std::map<WCoef, scalar_t> WilsonInterface::getAllFullRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    std::vector<WCoef> ids = WCoefMapper::get_group(group);
    std::map<WCoef, scalar_t> Cs;
    for (auto c : ids) {
        Cs.emplace(c, getFullRunCoefficient(group, c, order, cont_type));
    }

    return Cs;
}

std::map<WCoef, scalar_t> WilsonInterface::getAFR(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis) {
    return getAllFullRunCoefficients(group, order, cont_type);
}