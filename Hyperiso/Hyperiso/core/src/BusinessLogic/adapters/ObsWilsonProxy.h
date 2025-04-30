#ifndef OBSWILSONPROXY_H
#define OBSWILSONPROXY_H

#include "Include.h"
#include "IObsWilsonProxy.h"
#include "WilsonProvider.h"
#include "WilsonBuilder.h"

class ObsWilsonBuilder;

class ObsWilsonProxy : public IObsWilsonProxy<ObsWilsonBuilder> {
public:
    ObsWilsonProxy(std::shared_ptr<WilsonProvider> wilson_provider) : wil_p(wilson_provider) {}

    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override;
    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override;
    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override;
    complex_t getFR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override;
    std::map<QCDOrder, complex_t> getSM(WGroup group, WCoef coeff, bool sm_only=false);
    std::map<QCDOrder, complex_t> getSR(WGroup group, WCoef coeff, bool sm_only=false);
    std::map<WCoef, complex_t> getAM(WGroup group, QCDOrder order, bool sm_only=false);
    std::map<WCoef, complex_t> getAR(WGroup group, QCDOrder order, bool sm_only=false);
    std::map<WCoef, complex_t> getAFM(WGroup group, QCDOrder order, bool sm_only=false);
    std::map<WCoef, complex_t> getAFR(WGroup group, QCDOrder order, bool sm_only=false);

    std::shared_ptr<ObsWilsonBuilder> get_builder() override;

private:
    std::shared_ptr<WilsonProvider> wil_p;
};

#endif // __OBSWILSONPROXY_H__
