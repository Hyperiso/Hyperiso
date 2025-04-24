#ifndef OBSWILSONPROXY_H
#define OBSWILSONPROXY_H

#include "Include.h"
#include "IObsWilsonProxy.h"
#include "WilsonInterface.h"
#include "WilsonAdapter.h"

class ObsWilsonProxy : public IObsWilsonProxy {
public:
    ObsWilsonProxy() {
        if (!WilsonAdapter::built) {
            LOG_ERROR("InitializationError", "WilsonAdapter not built. Please call WilsonAdapter::build() before using ObsWilsonProxy.");
        }

        this->wi = wi;
    }

    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override {
        return this->wi.getM(group, coeff, order, sm_only);
    }

    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override {
        return this->wi.getFM(group, coeff, order, sm_only);
    }

    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override {
        return this->wi.getR(group, coeff, order, sm_only);
    }

    complex_t getFR(WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) override {
        return this->wi.getFR(group, coeff, order, sm_only);
    }

    std::map<QCDOrder, complex_t> getSM(WGroup group, WCoef coeff, bool sm_only=false) override {
        return this->wi.getSM(group, coeff, sm_only);
    }

    std::map<QCDOrder, complex_t> getSR(WGroup group, WCoef coeff, bool sm_only=false) override {
        return this->wi.getSR(group, coeff, sm_only);
    }

    std::map<WCoef, complex_t> getAM(WGroup group, QCDOrder order, bool sm_only=false) override {
        return this->wi.getAM(group, order, sm_only);
    }

    std::map<WCoef, complex_t> getAR(WGroup group, QCDOrder order, bool sm_only=false) override {
        return this->wi.getAR(group, order, sm_only);
    }

    std::map<WCoef, complex_t> getAFM(WGroup group, QCDOrder order, bool sm_only=false) override {
        return this->wi.getAFM(group, order, sm_only);
    }

    std::map<WCoef, complex_t> getAFR(WGroup group, QCDOrder order, bool sm_only=false) override {
        return this->wi.getAFR(group, order, sm_only);
    }

private:
    WilsonInterface wi;
};

#endif // __OBSWILSONPROXY_H__
