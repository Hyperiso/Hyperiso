#ifndef OBSWILSONPROXY_H
#define OBSWILSONPROXY_H

#include "Include.h"
#include "IObsWilsonProxy.h"
#include "WilsonProvider.h"
#include "WilsonBuilder.h"

class ObsWilsonBuilder;

class ObsWilsonProxy : public IObsWilsonProxy {
public:
    ObsWilsonProxy(std::shared_ptr<WilsonProvider> wilson_provider) : wil_p(wilson_provider) {}

    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;
    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;
    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;
    complex_t getFR(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;
    std::map<QCDOrder, complex_t> getSM(WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL);
    std::map<QCDOrder, complex_t> getSR(WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL);
    std::map<WCoef, complex_t> getAM(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);
    std::map<WCoef, complex_t> getAR(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);
    std::map<WCoef, complex_t> getAFM(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);
    std::map<WCoef, complex_t> getAFR(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);

    std::shared_ptr<IObsWilsonBuilder> get_builder() override;
    std::unordered_set<WilsonBasis> get_bases(WGroupId) override;
    void set_basis(WilsonBasis basis) override;

private:
    std::shared_ptr<WilsonProvider> wil_p;
    WilsonBasis basis {WilsonBasis::B_STANDARD}; 
};

#endif // __OBSWILSONPROXY_H__
