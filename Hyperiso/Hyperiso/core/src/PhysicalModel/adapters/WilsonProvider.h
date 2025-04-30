#ifndef __WILSONPROVIDER_H__
#define __WILSONPROVIDER_H__

#include "IWilsonProvider.h"
#include "Include.h"
#include "WilsonManager.h"

struct WilsonRequest : public AbstractConfig {
    WGroup group;
    WCoef coefficient;
    QCDOrder order {QCDOrder::LO};
    ContributionType contribution {ContributionType::TOTAL};
    ScaleType scale_type {ScaleType::HADRONIC};
    bool sum_qcd_orders {false};
};

class WilsonProvider : public IWilsonProvider {
public:
    WilsonProvider(std::shared_ptr<CoefficientManager> manager);

    scalar_t get(std::shared_ptr<AbstractConfig> config) override;

private:
    std::shared_ptr<CoefficientManager> cm;
};

#endif // __WILSONPROVIDER_H__
