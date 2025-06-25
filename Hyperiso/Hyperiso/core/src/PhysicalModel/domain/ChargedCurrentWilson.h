#ifndef __CHARGEDCURRENTWILSON_SUPER_H__
#define __CHARGEDCURRENTWILSON_SUPER_H__

#include "Wilson.h"

class C_V1 : public WilsonCoefficient {
public:
    C_V1();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&);
};

class C_V2 : public WilsonCoefficient {
public:
    C_V2() : WilsonCoefficient("C_V2", GroupMapper::str(WGroup::BCC, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2>(*this);
    }
};

class C_S1 : public WilsonCoefficient {
public:
    C_S1() : WilsonCoefficient("C_S1", GroupMapper::str(WGroup::BCC, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1>(*this);
    }

};  

class C_S2 : public WilsonCoefficient {
public:
    C_S2() : WilsonCoefficient("C_S2", GroupMapper::str(WGroup::BCC, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2>(*this);
    }
};

class C_T : public WilsonCoefficient {
public:
    C_T() : WilsonCoefficient("C_T", GroupMapper::str(WGroup::BCC, ScaleType::MATCHING)) {}

    
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T>(*this);
    }

};

#endif // __CHARGEDCURRENTWILSON_H__
