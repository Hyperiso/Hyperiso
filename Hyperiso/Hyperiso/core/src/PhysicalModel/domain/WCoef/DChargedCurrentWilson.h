#ifndef DCHARGEDCURRENTWILSON_H
#define DCHARGEDCURRENTWILSON_H

#include "Wilson.h"

class C_V1_cs : public WilsonCoefficient {
public:
    C_V1_cs();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_cs>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&);
};

class C_V2_cs : public WilsonCoefficient {
public:
    C_V2_cs() : WilsonCoefficient("C_V2_cs", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_cs>(*this);
    }
};

class C_S1_cs : public WilsonCoefficient {
public:
    C_S1_cs() : WilsonCoefficient("C_S1_cs", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_cs>(*this);
    }

};  

class C_S2_cs : public WilsonCoefficient {
public:
    C_S2_cs() : WilsonCoefficient("C_S2_cs", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_cs>(*this);
    }
};

class C_T_cs : public WilsonCoefficient {
public:
    C_T_cs() : WilsonCoefficient("C_T_cs", GroupMapper::str(WGroup::CC_cs, ScaleType::MATCHING)) {}

    
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_cs>(*this);
    }

};



class C_V1_cd : public WilsonCoefficient {
public:
    C_V1_cd();

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_cd>(*this);
    }

    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>&);
};

class C_V2_cd : public WilsonCoefficient {
public:
    C_V2_cd() : WilsonCoefficient("C_V2_cd", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_cd>(*this);
    }
};

class C_S1_cd : public WilsonCoefficient {
public:
    C_S1_cd() : WilsonCoefficient("C_S1_cd", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_cd>(*this);
    }

};  

class C_S2_cd : public WilsonCoefficient {
public:
    C_S2_cd() : WilsonCoefficient("C_S2_cd", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING)) {}

    
     
     

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_cd>(*this);
    }
};

class C_T_cd : public WilsonCoefficient {
public:
    C_T_cd() : WilsonCoefficient("C_T_cd", GroupMapper::str(WGroup::CC_cd, ScaleType::MATCHING)) {}

    
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_cd>(*this);
    }

};

#endif