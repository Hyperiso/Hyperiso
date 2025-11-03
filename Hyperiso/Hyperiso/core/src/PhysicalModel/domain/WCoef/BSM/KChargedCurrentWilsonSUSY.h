#ifndef KCHARGEDCURRENTWILSON_SUSY_H
#define KCHARGEDCURRENTWILSON_SUSY_H

#include "Wilson.h"

class C_V1_su_SUSY : public WilsonCoefficient {
public:
    C_V1_su_SUSY() : WilsonCoefficient("C_V1_su_SUSY", GroupMapper::str(WGroup::CC_su) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_su_SUSY>(*this);
    }

};

class C_V2_su_SUSY : public WilsonCoefficient {
public:
    C_V2_su_SUSY() : WilsonCoefficient("C_V2_su_SUSY", GroupMapper::str(WGroup::CC_su) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_su_SUSY>(*this);
    }

};

class C_S1_su_SUSY : public WilsonCoefficient {
public:
    C_S1_su_SUSY();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_su_SUSY>(*this);
    }

};

class C_S2_su_SUSY : public WilsonCoefficient {
public:
    C_S2_su_SUSY();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_su_SUSY>(*this);
    }

};

class C_T_su_SUSY : public WilsonCoefficient {
public:
    C_T_su_SUSY() : WilsonCoefficient("C_T_su_SUSY", GroupMapper::str(WGroup::CC_su) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_su_SUSY>(*this);
    }

};

#endif