#ifndef __CHARGEDCURRENTWILSON_SUPER_H__
#define __CHARGEDCURRENTWILSON_SUPER_H__

#include "WilsonSuper.h"

// class C_Blnu_A : public WilsonCoefficient {
// public:
//     C_Blnu_A() : WilsonCoefficient("C_Blnu_A", GroupMapper::str(WGroup::Blnu) + "_MATCH") {}

    
     
    

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C_Blnu_A>(*this);
//     }

//     static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&);
// };

// class C_Blnu_P : public WilsonCoefficient {
// public:
//     C_Blnu_P() : WilsonCoefficient("C_Blnu_P", GroupMapper::str(WGroup::Blnu) + "_MATCH") {}

    
     
     

//     std::shared_ptr<WilsonCoefficient> clone() const override {
//         return std::make_shared<C_Blnu_P>(*this);
//     }
// };

class C_V1 : public WilsonCoefficient {
public:
    C_V1() : WilsonCoefficient("C_V1", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&);
};

class C_V2 : public WilsonCoefficient {
public:
    C_V2() : WilsonCoefficient("C_V2", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2>(*this);
    }
};

class C_S1 : public WilsonCoefficient {
public:
    C_S1() : WilsonCoefficient("C_S1", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1>(*this);
    }

};  

class C_S2 : public WilsonCoefficient {
public:
    C_S2() : WilsonCoefficient("C_S2", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2>(*this);
    }
};

class C_T : public WilsonCoefficient {
public:
    C_T() : WilsonCoefficient("C_T", GroupMapper::str(WGroup::BCC) + "_MATCH") {}

    
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T>(*this);
    }

};

#endif // __CHARGEDCURRENTWILSON_H__
