#ifndef KCHARGEDCURRENTWILSON_H
#define KCHARGEDCURRENTWILSON_H

#include "Wilson.h"

class C_V1_su : public WilsonCoefficient {
public:
    C_V1_su();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_su>(*this);
    }

    static double compute_LO(const ParamSrc&);
};

class C_V2_su : public WilsonCoefficient {
public:
    C_V2_su() : WilsonCoefficient("C_V2_su", GroupMapper::str(WGroup::CC_su, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_su>(*this);
    }
};

class C_S1_su : public WilsonCoefficient {
public:
    C_S1_su() : WilsonCoefficient("C_S1_su", GroupMapper::str(WGroup::CC_su, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_su>(*this);
    }

};  

class C_S2_su : public WilsonCoefficient {
public:
    C_S2_su() : WilsonCoefficient("C_S2_su", GroupMapper::str(WGroup::CC_su, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_su>(*this);
    }
};

class C_T_su : public WilsonCoefficient {
public:
    C_T_su() : WilsonCoefficient("C_T_su", GroupMapper::str(WGroup::CC_su, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_su>(*this);
    }
};

#endif