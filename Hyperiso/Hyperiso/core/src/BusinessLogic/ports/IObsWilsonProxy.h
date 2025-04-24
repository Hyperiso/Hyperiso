#ifndef __IOBSWILSONPROXY_H__
#define __IOBSWILSONPROXY_H__

#include "Include.h"

class IObsWilsonProxy {
public:
    virtual ~IObsWilsonProxy() = default;

    virtual complex_t                       getM    (WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) = 0;
    virtual complex_t                       getFM   (WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) = 0;
    virtual complex_t                       getR    (WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) = 0;
    virtual complex_t                       getFR   (WGroup group, WCoef coeff, QCDOrder order, bool sm_only=false) = 0;
    virtual std::map<QCDOrder, complex_t>   getSM   (WGroup group, WCoef coeff, bool sm_only=false)                 = 0;
    virtual std::map<QCDOrder, complex_t>   getSR   (WGroup group, WCoef coeff, bool sm_only=false)                 = 0;
    virtual std::map<WCoef, complex_t>      getAM   (WGroup group, QCDOrder order, bool sm_only=false)              = 0;
    virtual std::map<WCoef, complex_t>      getAR   (WGroup group, QCDOrder order, bool sm_only=false)              = 0;
    virtual std::map<WCoef, complex_t>      getAFM  (WGroup group, QCDOrder order, bool sm_only=false)              = 0;
    virtual std::map<WCoef, complex_t>      getAFR  (WGroup group, QCDOrder order, bool sm_only=false)              = 0;
};

#endif // __IOBSWILSONPROXY_H__
