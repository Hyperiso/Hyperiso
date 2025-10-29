#ifndef CHARGEDCURRENTWILSON_H
#define CHARGEDCURRENTWILSON_H

#include "Wilson.h"

class C_V1_bc : public WilsonCoefficient {
public:
    C_V1_bc();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_bc>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&);
};

class C_V2_bc : public WilsonCoefficient {
public:
    C_V2_bc() : WilsonCoefficient("C_V2_bc", GroupMapper::str(WGroup::BCC_bc, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_bc>(*this);
    }
};

class C_S1_bc : public WilsonCoefficient {
public:
    C_S1_bc() : WilsonCoefficient("C_S1_bc", GroupMapper::str(WGroup::BCC_bc, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_bc>(*this);
    }

};  

class C_S2_bc : public WilsonCoefficient {
public:
    C_S2_bc() : WilsonCoefficient("C_S2_bc", GroupMapper::str(WGroup::BCC_bc, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_bc>(*this);
    }
};

class C_T_bc : public WilsonCoefficient {
public:
    C_T_bc() : WilsonCoefficient("C_T_bc", GroupMapper::str(WGroup::BCC_bc, ScaleType::MATCHING)) {}

    
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_bc>(*this);
    }

};



class C_V1_bu : public WilsonCoefficient {
public:
    C_V1_bu();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_bu>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&);
};

class C_V2_bu : public WilsonCoefficient {
public:
    C_V2_bu() : WilsonCoefficient("C_V2_bu", GroupMapper::str(WGroup::BCC_bu, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_bu>(*this);
    }
};

class C_S1_bu : public WilsonCoefficient {
public:
    C_S1_bu() : WilsonCoefficient("C_S1_bu", GroupMapper::str(WGroup::BCC_bu, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_bu>(*this);
    }

};  

class C_S2_bu : public WilsonCoefficient {
public:
    C_S2_bu() : WilsonCoefficient("C_S2_bu", GroupMapper::str(WGroup::BCC_bu, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_bu>(*this);
    }
};

class C_T_bu : public WilsonCoefficient {
public:
    C_T_bu() : WilsonCoefficient("C_T_bu", GroupMapper::str(WGroup::BCC_bu, ScaleType::MATCHING)) {}

    
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_bu>(*this);
    }

};

#endif // __CHARGEDCURRENTWILSON_H__
