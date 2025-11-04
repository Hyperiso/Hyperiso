#ifndef BSCALAR_WILSON_H
#define BSCALAR_WILSON_H

#include "Wilson.h"

class CQ1 : public WilsonCoefficient {
public:
    CQ1();
    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ1>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CQ2 : public WilsonCoefficient {
public:
    CQ2();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ2>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

#endif