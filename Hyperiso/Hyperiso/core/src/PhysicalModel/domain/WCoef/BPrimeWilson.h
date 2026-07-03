#ifndef BPRIME_WILSON_H
#define BPRIME_WILSON_H

#include "Wilson.h"
#include "wcoef_ids.hpp"

class CP1 : public WilsonCoefficient {
public:
    CP1() : WilsonCoefficient("CP1", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP1>(*this);
    }

};

class CP2 : public WilsonCoefficient {
public:
    CP2() : WilsonCoefficient("CP2", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP2>(*this);
    }

};

class CP3 : public WilsonCoefficient {
public:
    CP3() : WilsonCoefficient("CP3", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP3>(*this);
    }

};

class CP4 : public WilsonCoefficient {
public:
    CP4() : WilsonCoefficient("CP4", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP4>(*this);
    }

};

class CP5 : public WilsonCoefficient {
public:
    CP5() : WilsonCoefficient("CP5", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP5>(*this);
    }

};

class CP6 : public WilsonCoefficient {
public:
    CP6() : WilsonCoefficient("CP6", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP6>(*this);
    }

};

class CP7 : public WilsonCoefficient {
public:
    CP7();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP7>(*this);
    }
    static double compute_LO(const ParamSrc& src);
};

class CP8 : public WilsonCoefficient {
public:
    CP8();

    
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP8>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class CP9 : public WilsonCoefficient {
public:
    CP9() : WilsonCoefficient("CP9", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP9>(*this);
    }

};

class CP10 : public WilsonCoefficient {
public:
    CP10() : WilsonCoefficient("CP10", GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CP10>(*this);
    }

};

class CPQ1 : public WilsonCoefficient {
public:
    explicit CPQ1(WCoef coef = WCoef::CPQ1_MU) : WilsonCoefficient(WCoefMapper::str(coef), GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ1>(*this);
    }
};

class CPQ2 : public WilsonCoefficient {
public:
    explicit CPQ2(WCoef coef = WCoef::CPQ2_MU) : WilsonCoefficient(WCoefMapper::str(coef), GroupMapper::str(WGroup::BPrime) + "_MATCH") {}

     
     
    void NNLO_calculation() {} 

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CPQ2>(*this);
    }
};


#endif