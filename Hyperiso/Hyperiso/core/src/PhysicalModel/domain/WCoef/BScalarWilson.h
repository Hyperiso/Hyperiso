#ifndef BSCALAR_WILSON_H
#define BSCALAR_WILSON_H

#include "Wilson.h"
#include "wcoef_ids.hpp"

class CQ1 : public WilsonCoefficient {
public:
    explicit CQ1(WCoef coef = WCoef::CQ1_MU);
    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ1>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_LO(const ParamSrc& src, int lepton_mass_slot);
};

class CQ2 : public WilsonCoefficient {
public:
    explicit CQ2(WCoef coef = WCoef::CQ2_MU);

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CQ2>(*this);
    }

    static double compute_LO(const ParamSrc& src);
    static double compute_LO(const ParamSrc& src, int lepton_mass_slot);
};

#endif
