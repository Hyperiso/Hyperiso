#ifndef PICHARGEDCURRENTWILSON_H
#define PICHARGEDCURRENTWILSON_H

#include "Wilson.h"

class C_V1_du : public WilsonCoefficient {
public:
    C_V1_du();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_du>(*this);
    }

    static double compute_LO(const ParamSrc& src);
};

class C_V2_du : public WilsonCoefficient {
public:
    C_V2_du() : WilsonCoefficient("C_V2_du", GroupMapper::str(WGroup::CC_du, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_du>(*this);
    }
};

class C_S1_du : public WilsonCoefficient {
public:
    C_S1_du() : WilsonCoefficient("C_S1_du", GroupMapper::str(WGroup::CC_du, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_du>(*this);
    }

};  

class C_S2_du : public WilsonCoefficient {
public:
    C_S2_du() : WilsonCoefficient("C_S2_du", GroupMapper::str(WGroup::CC_du, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_du>(*this);
    }
};

class C_T_du : public WilsonCoefficient {
public:
    C_T_du() : WilsonCoefficient("C_T_du", GroupMapper::str(WGroup::CC_du, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_du>(*this);
    }
};

#endif