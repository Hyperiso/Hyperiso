#include "ObsWilsonProxy.h"
#include "ObsWilsonBuilder.h"

complex_t ObsWilsonProxy::getM(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution) {
    WilsonRequest request{group, coeff, order, contribution, ScaleType::MATCHING, false};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getFM(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution) {
    WilsonRequest request{group, coeff, order, contribution, ScaleType::MATCHING, true};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getR(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution) {
    WilsonRequest request{group, coeff, order, contribution, ScaleType::HADRONIC, false};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getFR(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution){
    WilsonRequest request{group, coeff, order, contribution, ScaleType::HADRONIC, true};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getM(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution) {
    WilsonRequest request{group, coeff, order, contribution, ScaleType::MATCHING, false};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getFM(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution) {
    WilsonRequest request{group, coeff, order, contribution, ScaleType::MATCHING, true};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getR(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution) {
    WilsonRequest request{group, coeff, order, contribution, ScaleType::HADRONIC, false};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getFR(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType contribution){
    WilsonRequest request{group, coeff, order, contribution, ScaleType::HADRONIC, true};
    request.basis = this->basis;
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

std::map<QCDOrder, complex_t> ObsWilsonProxy::getSM(WGroup group, WCoef coeff, ContributionType contribution) {
    return {
        {QCDOrder::LO, getM(group, coeff, QCDOrder::LO, contribution)},
        {QCDOrder::NLO, getM(group, coeff, QCDOrder::NLO, contribution)},
        {QCDOrder::NNLO, getM(group, coeff, QCDOrder::NNLO, contribution)}
    };
}

std::map<QCDOrder, complex_t> ObsWilsonProxy::getSR(WGroup group, WCoef coeff, ContributionType contribution) {
    return {
        {QCDOrder::LO, getR(group, coeff, QCDOrder::LO, contribution)},
        {QCDOrder::NLO, getR(group, coeff, QCDOrder::NLO, contribution)},
        {QCDOrder::NNLO, getR(group, coeff, QCDOrder::NNLO, contribution)}
    };
}

std::map<QCDOrder, complex_t> ObsWilsonProxy::getSM(WGroupId group, WCoefId coeff, ContributionType contribution) {
    return {
        {QCDOrder::LO, getM(group, coeff, QCDOrder::LO, contribution)},
        {QCDOrder::NLO, getM(group, coeff, QCDOrder::NLO, contribution)},
        {QCDOrder::NNLO, getM(group, coeff, QCDOrder::NNLO, contribution)}
    };
}

std::map<QCDOrder, complex_t> ObsWilsonProxy::getSR(WGroupId group, WCoefId coeff, ContributionType contribution) {
    return {
        {QCDOrder::LO, getR(group, coeff, QCDOrder::LO, contribution)},
        {QCDOrder::NLO, getR(group, coeff, QCDOrder::NLO, contribution)},
        {QCDOrder::NNLO, getR(group, coeff, QCDOrder::NNLO, contribution)}
    };
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAM(WGroup group, QCDOrder order, ContributionType contribution) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getM(group, Ci, order, contribution);
    }
    return Cs;
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAR(WGroup group, QCDOrder order, ContributionType contribution) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getR(group, Ci, order, contribution);
    }
    return Cs;
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAFM(WGroup group, QCDOrder order, ContributionType contribution) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getFM(group, Ci, order, contribution);
    }
    return Cs;
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAFR(WGroup group, QCDOrder order, ContributionType contribution) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getFR(group, Ci, order, contribution);
    }
    return Cs;
}

std::shared_ptr<IObsWilsonBuilder> ObsWilsonProxy::get_builder() {
    auto wilson_builder = std::static_pointer_cast<WilsonBuilder>(this->wil_p->get_builder());
    return std::make_shared<ObsWilsonBuilder>(wilson_builder);
}

std::unordered_set<WilsonBasis> ObsWilsonProxy::get_bases(WGroupId group) {
    return this->wil_p->get_bases(group);
}

void ObsWilsonProxy::set_basis(WilsonBasis basis) {
    this->basis = basis;
}
