#ifndef DCHARGEDCURRENTWILSON_SUSY_H
#define DCHARGEDCURRENTWILSON_SUSY_H

#include "Wilson.h"

class C_V1_cs_SUSY : public WilsonCoefficient {
public:
    C_V1_cs_SUSY() : WilsonCoefficient("C_V1_cs_SUSY", GroupMapper::str(WGroup::CC_cs) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_cs_SUSY>(*this);
    }

};

class C_V2_cs_SUSY : public WilsonCoefficient {
public:
    C_V2_cs_SUSY() : WilsonCoefficient("C_V2_cs_SUSY", GroupMapper::str(WGroup::CC_cs) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_cs_SUSY>(*this);
    }

};

class C_S1_cs_SUSY : public WilsonCoefficient {
public:
    C_S1_cs_SUSY();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_cs_SUSY>(*this);
    }

};

class C_S2_cs_SUSY : public WilsonCoefficient {
public:
    C_S2_cs_SUSY();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_cs_SUSY>(*this);
    }

};

class C_T_cs_SUSY : public WilsonCoefficient {
public:
    C_T_cs_SUSY() : WilsonCoefficient("C_T_cs_SUSY", GroupMapper::str(WGroup::CC_cs) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_cs_SUSY>(*this);
    }

};



class C_V1_cd_SUSY : public WilsonCoefficient {
public:
    C_V1_cd_SUSY() : WilsonCoefficient("C_V1_cd_SUSY", GroupMapper::str(WGroup::CC_cd) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_cd_SUSY>(*this);
    }

};

class C_V2_cd_SUSY : public WilsonCoefficient {
public:
    C_V2_cd_SUSY() : WilsonCoefficient("C_V2_cd_SUSY", GroupMapper::str(WGroup::CC_cd) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_cd_SUSY>(*this);
    }

};

class C_S1_cd_SUSY : public WilsonCoefficient {
public:
    C_S1_cd_SUSY();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_cd_SUSY>(*this);
    }

};

class C_S2_cd_SUSY : public WilsonCoefficient {
public:
    C_S2_cd_SUSY();

    
     
    static double compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_cd_SUSY>(*this);
    }

};

class C_T_cd_SUSY : public WilsonCoefficient {
public:
    C_T_cd_SUSY() : WilsonCoefficient("C_T_cd_SUSY", GroupMapper::str(WGroup::CC_cd) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_cd_SUSY>(*this);
    }

};

#endif