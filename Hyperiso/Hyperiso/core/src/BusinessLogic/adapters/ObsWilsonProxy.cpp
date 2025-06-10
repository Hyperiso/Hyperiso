#include "ObsWilsonProxy.h"
#include "ObsWilsonBuilder.h"

complex_t ObsWilsonProxy::getM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only) {
    WilsonRequest request{group, coeff, order, sm_only ? ContributionType::SM : ContributionType::TOTAL, ScaleType::MATCHING, false};
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getFM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only) {
    WilsonRequest request{group, coeff, order, sm_only ? ContributionType::SM : ContributionType::TOTAL, ScaleType::MATCHING, true};
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only) {
    WilsonRequest request{group, coeff, order, sm_only ? ContributionType::SM : ContributionType::TOTAL, ScaleType::HADRONIC, false};
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

complex_t ObsWilsonProxy::getFR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only){
    WilsonRequest request{group, coeff, order, sm_only ? ContributionType::SM : ContributionType::TOTAL, ScaleType::HADRONIC, true};
    return this->wil_p->get(std::make_shared<WilsonRequest>(request));
}

std::map<QCDOrder, complex_t> ObsWilsonProxy::getSM(WGroup group, WCoef coeff, bool sm_only) {
    return {
        {QCDOrder::LO, getM(group, coeff, QCDOrder::LO, sm_only)},
        {QCDOrder::NLO, getM(group, coeff, QCDOrder::NLO, sm_only)},
        {QCDOrder::NNLO, getM(group, coeff, QCDOrder::NNLO, sm_only)}
    };
}

std::map<QCDOrder, complex_t> ObsWilsonProxy::getSR(WGroup group, WCoef coeff, bool sm_only) {
    return {
        {QCDOrder::LO, getR(group, coeff, QCDOrder::LO, sm_only)},
        {QCDOrder::NLO, getR(group, coeff, QCDOrder::NLO, sm_only)},
        {QCDOrder::NNLO, getR(group, coeff, QCDOrder::NNLO, sm_only)}
    };
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAM(WGroup group, QCDOrder order, bool sm_only) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getM(group, Ci, order, sm_only);
    }
    return Cs;
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAR(WGroup group, QCDOrder order, bool sm_only) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getR(group, Ci, order, sm_only);
    }
    return Cs;
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAFM(WGroup group, QCDOrder order, bool sm_only) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getFM(group, Ci, order, sm_only);
    }
    return Cs;
}

std::map<WCoef, complex_t> ObsWilsonProxy::getAFR(WGroup group, QCDOrder order, bool sm_only) {
    std::map<WCoef, complex_t> Cs;
    for (WCoef Ci : WCoefMapper::get_group(group)) {
        Cs[Ci] = getFR(group, Ci, order, sm_only);
    }
    return Cs;
}

std::shared_ptr<ObsWilsonBuilder> ObsWilsonProxy::get_builder() {
    auto wilson_builder = std::static_pointer_cast<WilsonBuilder>(this->wil_p->get_builder());
    return std::make_shared<ObsWilsonBuilder>(wilson_builder);
}

std::unordered_set<WilsonBasis> ObsWilsonProxy::get_bases(WGroup group) {
    return this->wil_p->get_bases(group);
}