#ifndef PICHARGEDCURRENTWILSON_SUSY_H
#define PICHARGEDCURRENTWILSON_SUSY_H

#include "Wilson.h"

class C_V1_du_SUSY : public WilsonCoefficient {
public:
    C_V1_du_SUSY() : WilsonCoefficient("C_V1_du_SUSY", GroupMapper::str(WGroup::CC_du) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_du_SUSY>(*this);
    }

};

class C_V2_du_SUSY : public WilsonCoefficient {
public:
    C_V2_du_SUSY() : WilsonCoefficient("C_V2_du_SUSY", GroupMapper::str(WGroup::CC_du) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_du_SUSY>(*this);
    }

};

class C_S1_du_SUSY : public WilsonCoefficient {
public:
    C_S1_du_SUSY();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_du_SUSY>(*this);
    }

};

class C_S2_du_SUSY : public WilsonCoefficient {
public:
    C_S2_du_SUSY();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_du_SUSY>(*this);
    }

};

class C_T_du_SUSY : public WilsonCoefficient {
public:
    C_T_du_SUSY() : WilsonCoefficient("C_T_du_SUSY", GroupMapper::str(WGroup::CC_du) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_du_SUSY>(*this);
    }

};

#endif