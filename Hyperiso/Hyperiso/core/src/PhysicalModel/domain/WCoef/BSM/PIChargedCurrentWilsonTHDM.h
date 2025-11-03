#ifndef PICHARGEDCURRENTWILSON_THDM_H
#define PIHARGEDCURRENTWILSON_THDM_H

#include "Wilson.h"

class C_V1_du_THDM : public WilsonCoefficient {
public:
    C_V1_du_THDM() : WilsonCoefficient("C_V1_du_THDM", GroupMapper::str(WGroup::CC_du) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_du_THDM>(*this);
    }

};

class C_V2_du_THDM : public WilsonCoefficient {
public:
    C_V2_du_THDM() : WilsonCoefficient("C_V2_du_THDM", GroupMapper::str(WGroup::CC_du) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_du_THDM>(*this);
    }

};

class C_S1_du_THDM : public WilsonCoefficient {
public:
    C_S1_du_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_du_THDM>(*this);
    }

};

class C_S2_du_THDM : public WilsonCoefficient {
public:
    C_S2_du_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_du_THDM>(*this);
    }

};

class C_T_du_THDM : public WilsonCoefficient {
public:
    C_T_du_THDM() : WilsonCoefficient("C_T_du_THDM", GroupMapper::str(WGroup::CC_du) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_du_THDM>(*this);
    }

};

#endif