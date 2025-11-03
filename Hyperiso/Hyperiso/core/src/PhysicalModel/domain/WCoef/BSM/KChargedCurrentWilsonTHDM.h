#ifndef KCHARGEDCURRENTWILSON_THDM_H
#define KCHARGEDCURRENTWILSON_THDM_H

#include "Wilson.h"

class C_V1_su_THDM : public WilsonCoefficient {
public:
    C_V1_su_THDM() : WilsonCoefficient("C_V1_su_THDM", GroupMapper::str(WGroup::CC_su) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_su_THDM>(*this);
    }

};

class C_V2_su_THDM : public WilsonCoefficient {
public:
    C_V2_su_THDM() : WilsonCoefficient("C_V2_su_THDM", GroupMapper::str(WGroup::CC_su) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_su_THDM>(*this);
    }

};

class C_S1_su_THDM : public WilsonCoefficient {
public:
    C_S1_su_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_su_THDM>(*this);
    }

};

class C_S2_su_THDM : public WilsonCoefficient {
public:
    C_S2_su_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_su_THDM>(*this);
    }

};

class C_T_su_THDM : public WilsonCoefficient {
public:
    C_T_su_THDM() : WilsonCoefficient("C_T_su_THDM", GroupMapper::str(WGroup::CC_su) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_su_THDM>(*this);
    }

};

#endif