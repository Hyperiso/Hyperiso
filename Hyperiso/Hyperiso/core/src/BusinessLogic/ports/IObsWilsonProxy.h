#ifndef IOBSWILSONPROXY_H
#define IOBSWILSONPROXY_H

#include "Include.h"

template<typename BuilderType>
class IObsWilsonProxy {
public:
    virtual ~IObsWilsonProxy() = default;

    virtual complex_t                       getM    (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;
    virtual complex_t                       getFM   (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;
    virtual complex_t                       getR    (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;
    virtual complex_t                       getFR   (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;
    virtual std::map<QCDOrder, complex_t>   getSM   (WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL)                 = 0;
    virtual std::map<QCDOrder, complex_t>   getSR   (WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL)                 = 0;
    virtual std::map<WCoef, complex_t>      getAM   (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;
    virtual std::map<WCoef, complex_t>      getAR   (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;
    virtual std::map<WCoef, complex_t>      getAFM  (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;
    virtual std::map<WCoef, complex_t>      getAFR  (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;

    virtual std::shared_ptr<BuilderType> get_builder() = 0;
    virtual std::unordered_set<WilsonBasis> get_bases(WGroupId) = 0;
};

#endif // __IOBSWILSONPROXY_H__
