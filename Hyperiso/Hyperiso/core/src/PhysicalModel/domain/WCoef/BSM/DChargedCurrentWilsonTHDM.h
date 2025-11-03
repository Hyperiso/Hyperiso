#ifndef DCHARGEDCURRENTWILSON_THDM_H
#define DCHARGEDCURRENTWILSON_THDM_H

#include "Wilson.h"

class C_V1_cs_THDM : public WilsonCoefficient {
public:
    C_V1_cs_THDM() : WilsonCoefficient("C_V1_cs_THDM", GroupMapper::str(WGroup::CC_cs) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_cs_THDM>(*this);
    }

};

class C_V2_cs_THDM : public WilsonCoefficient {
public:
    C_V2_cs_THDM() : WilsonCoefficient("C_V2_cs_THDM", GroupMapper::str(WGroup::CC_cs) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_cs_THDM>(*this);
    }

};

class C_S1_cs_THDM : public WilsonCoefficient {
public:
    C_S1_cs_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_cs_THDM>(*this);
    }

};

class C_S2_cs_THDM : public WilsonCoefficient {
public:
    C_S2_cs_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_cs_THDM>(*this);
    }

};

class C_T_cs_THDM : public WilsonCoefficient {
public:
    C_T_cs_THDM() : WilsonCoefficient("C_T_cs_THDM", GroupMapper::str(WGroup::CC_cs) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_cs_THDM>(*this);
    }

};



class C_V1_cd_THDM : public WilsonCoefficient {
public:
    C_V1_cd_THDM() : WilsonCoefficient("C_V1_cd_THDM", GroupMapper::str(WGroup::CC_cd) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V1_cd_THDM>(*this);
    }

};

class C_V2_cd_THDM : public WilsonCoefficient {
public:
    C_V2_cd_THDM() : WilsonCoefficient("C_V2_cd_THDM", GroupMapper::str(WGroup::CC_cd) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_V2_cd_THDM>(*this);
    }

};

class C_S1_cd_THDM : public WilsonCoefficient {
public:
    C_S1_cd_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S1_cd_THDM>(*this);
    }

};

class C_S2_cd_THDM : public WilsonCoefficient {
public:
    C_S2_cd_THDM();

    
     
    static double compute_LO(const ParamSrc& src);

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_S2_cd_THDM>(*this);
    }

};

class C_T_cd_THDM : public WilsonCoefficient {
public:
    C_T_cd_THDM() : WilsonCoefficient("C_T_cd_THDM", GroupMapper::str(WGroup::CC_cd) + "_MATCH") {}

     
     
    

    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<C_T_cd_THDM>(*this);
    }

};

#endif